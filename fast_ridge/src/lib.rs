use core::f64;

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{Read, Write};

#[derive(Clone, Debug)]
struct RidgeResult {
    gcv_score: f64,
    best_alpha: f64,
    bias: f64,
    weights: DVector<f64>,
}

fn build_mat_with_edges(x_nodes: &DMatrix<f64>, selected_edges: &[(usize, usize)]) -> DMatrix<f64> {
    let (m, k) = x_nodes.shape();
    let l = selected_edges.len();
    let d = k + l;

    let mut x_final = DMatrix::zeros(m, d);

    //copy nodes features
    for c in 0..k {
        let col = x_nodes.column(c);
        x_final.column_mut(c).copy_from(&col);
    }

    // build edges features
    for (idx, &(u, v)) in selected_edges.iter().enumerate() {
        let col_idx = k + idx;
        for r in 0..m {
            if x_nodes[(r, u)] > 0.5 && x_nodes[(r, v)] > 0.5 {
                x_final[(r, col_idx)] = 1.0;
            }
        }
    }
    x_final
}
fn ridge_gcv(x: &DMatrix<f64>, y: &DVector<f64>, alphas: &[f64]) -> RidgeResult {
    let (m, _d) = x.shape();

    // Centering y
    let y_mean = y.mean();
    let yc = y.add_scalar(-y_mean);

    let mu = x.row_mean().transpose();
    let kernel_mat = x * x.transpose();

    let ones = DVector::from_element(m, 1.0);
    let v = x * &mu;
    let s = mu.dot(&mu);

    let v_one = &v * ones.transpose();
    let one_v = &ones * v.transpose();
    let one_one = &ones * ones.transpose();

    let kernel_centered = kernel_mat - v_one - one_v + one_one.scale(s);

    let eigen = SymmetricEigen::new(kernel_centered.clone());
    let lam = eigen.eigenvalues;
    let u_mat = eigen.eigenvectors;
    let uy = u_mat.transpose() * &yc;

    let mut best_gcv = f64::INFINITY;
    let mut best_alpha = alphas[0];

    for &alpha in alphas {
        let denom = lam.add_scalar(alpha);
        let r_vec: DVector<f64> = uy.zip_map(&denom, |uyi, den| (alpha / den) * uyi);
        let rss = r_vec.dot(&r_vec);
        let tr_h: f64 = lam.zip_map(&denom, |l, d| l / d).sum();
        let denom_gcv = (m as f64 - tr_h).max(1e-12).powi(2);
        let gcv = rss / denom_gcv;

        if gcv < best_gcv {
            best_gcv = gcv;
            best_alpha = alpha;
        }
    }
    let eye = DMatrix::identity(m, m);
    let lhs = kernel_centered + eye.scale(best_alpha);
    let a_vec = lhs.lu().solve(&yc).unwrap_or_else(|| DVector::zeros(m));

    let mut w = x.transpose() * &a_vec;
    let a_sum = a_vec.sum();
    w = w - mu.scale(a_sum);
    let b = y_mean - mu.dot(&w);

    RidgeResult {
        gcv_score: best_gcv,
        best_alpha,
        bias: b,
        weights: w,
    }
}

#[pyclass]
#[derive(Debug, Serialize, Deserialize)]
struct GraphModel {
    node_weights: Option<DVector<f64>>,
    edge_weights: Option<DVector<f64>>,
    selected_edges: Vec<(usize, usize)>,
    bias: f64,
    node_cnt: usize,
}

#[pymethods]
impl GraphModel {
    #[new]
    fn new() -> Self {
        GraphModel {
            node_weights: None,
            edge_weights: None,
            selected_edges: Vec::new(),
            bias: 0.0,
            node_cnt: 0,
        }
    }

    fn fit(
        &mut self,
        x_from_py: PyReadonlyArray2<'_, i8>,
        y_from_py: PyReadonlyArray1<'_, f64>,
        max_edges: usize,
        step_size: usize,
        min_edge_cnt: i32,
        threads: usize,
    ) -> PyResult<(f64, f64, usize)> {
        let x_np = x_from_py.as_array();
        let y_np = y_from_py.as_array();

        let (m, nc) = (x_np.shape()[0], x_np.shape()[1]);
        self.node_cnt = nc;

        let data_vec: Vec<f64> = x_np.iter().map(|&v| v as f64).collect();
        let x_nodes = DMatrix::from_row_slice(m, nc, &data_vec);
        let y = DVector::from_column_slice(y_np.as_slice().unwrap());

        let y_var = y.variance();
        let base_scale = y_var.max(1e-6);

        // generate log-uniform distribute in range of [-1000, 1000]*base_scale
        let alphas: Vec<f64> = (0..25)
            .map(|i| 10f64.powf(-3.0 + (6.0 * i as f64 / 24.0)) * base_scale)
            .collect();

        // Stage 1, modeling with nodes
        let res1 = ridge_gcv(&x_nodes, &y, &alphas);
        let stage1_alpha = res1.best_alpha;

        // Edge screening
        let preds = &x_nodes * &res1.weights;
        let resid = &y - preds.add_scalar(res1.bias);
        let resid_mean = resid.mean();
        let resid_centered = resid.add_scalar(-resid_mean);
        let resid_var = resid_centered.dot(&resid_centered) / (m as f64);

        let mut all_candidates = Vec::new();

        let num_threads = threads.max(1);
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("Failed to build thread pool");

        pool.install(|| {
            if resid_var > 1e-12 && max_edges > 0 {
                let num_pairs = nc * (nc - 1) / 2;

                let (final_counts, final_sums) = (0..m)
                    .into_par_iter()
                    .fold(
                        || (vec![0i32; num_pairs], vec![0.0f64; num_pairs]),
                        |(mut counts, mut sums), r| {
                            let mut active = Vec::with_capacity(nc);
                            for c in 0..nc {
                                if x_nodes[(r, c)] > 0.5 {
                                    active.push(c);
                                }
                            }
                            let rm = resid_centered[r];
                            let s = active.len();

                            if s >= 2 {
                                for i in 0..s - 1 {
                                    for j in i + 1..s {
                                        let u = active[i];
                                        let v = active[j];
                                        let (min_u, max_v) = if u < v { (u, v) } else { (v, u) };
                                        let idx =
                                            min_u * (2 * nc - min_u - 1) / 2 + (max_v - min_u - 1);
                                        counts[idx] += 1;
                                        sums[idx] += rm;
                                    }
                                }
                            }
                            (counts, sums)
                        },
                    )
                    .reduce(
                        || (vec![0; num_pairs], vec![0.0; num_pairs]),
                        |(mut c1, mut s1), (c2, s2)| {
                            for i in 0..num_pairs {
                                c1[i] += c2[i];
                                s1[i] += s2[i];
                            }
                            (c1, s1)
                        },
                    );

                for idx in 0..num_pairs {
                    if final_counts[idx] >= min_edge_cnt {
                        let p = final_counts[idx] as f64 / m as f64;
                        let z_var = p * (1.0 - p);
                        let cov = final_sums[idx] / m as f64;
                        let corr = cov / (z_var * resid_var).sqrt().max(1e-12);
                        all_candidates.push((idx, corr.abs()));
                    }
                }
                all_candidates.par_sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            }
        });

        // test edges [step, 2*step, ..., max_edges]
        let mut candidate_counts = Vec::new();
        if step_size > 0 {
            let mut curr = step_size;
            while curr <= max_edges && curr <= all_candidates.len() {
                candidate_counts.push(curr);
                curr += step_size;
            }
        }

        let mut best_res = res1;
        let mut best_cnt = 0;

        pool.install(|| {
            if !candidate_counts.is_empty() {
                let best_parallel_result = candidate_counts
                    .par_iter()
                    .map(|&k_cnt| {
                        let current_edges: Vec<(usize, usize)> = all_candidates
                            .iter()
                            .take(k_cnt)
                            .map(|&(idx, _)| {
                                let mut rem = idx;
                                for i in 0..nc - 1 {
                                    let blk = nc - 1 - i;
                                    if rem < blk {
                                        return (i, i + 1 + rem);
                                    }
                                    rem -= blk;
                                }
                                (0, 0)
                            })
                            .collect();
                        // build matrix with edges
                        let x_curr = build_mat_with_edges(&x_nodes, &current_edges);

                        // solve gcv
                        let res = ridge_gcv(&x_curr, &y, &alphas);

                        (res, k_cnt, current_edges)
                    })
                    .reduce(
                        || {
                            let dummy = RidgeResult {
                                gcv_score: f64::INFINITY,
                                best_alpha: 0.9,
                                bias: 0.0,
                                weights: DVector::zeros(0),
                            };
                            (dummy, 0, Vec::new())
                        },
                        |best_so_far, current| {
                            if current.0.gcv_score < best_so_far.0.gcv_score {
                                current
                            } else {
                                best_so_far
                            }
                        },
                    );
                if best_parallel_result.0.gcv_score < best_res.gcv_score {
                    best_res = best_parallel_result.0;
                    best_cnt = best_parallel_result.1;
                    self.selected_edges = best_parallel_result.2;
                }
            }
        });

        if best_cnt == 0 {
            self.selected_edges.clear();
        }

        self.bias = best_res.bias;
        self.node_weights = Some(best_res.weights.rows(0, nc).into_owned());
        if best_res.weights.len() > nc {
            self.edge_weights = Some(
                best_res
                    .weights
                    .rows(nc, best_res.weights.len() - nc)
                    .into_owned(),
            );
        } else {
            self.edge_weights = Some(DVector::zeros(0));
        }

        Ok((stage1_alpha, best_res.best_alpha, best_cnt))
    }

    fn predict<'py>(
        &self,
        py: Python<'py>,
        x_from_py: PyReadonlyArray2<'_, i8>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_np = x_from_py.as_array();
        let (m, k) = (x_np.shape()[0], x_np.shape()[1]);

        if k != self.node_cnt {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Input node count={} does not match Model node count={}",
                k, self.node_cnt
            )));
        }

        let data_vec: Vec<f64> = x_np.iter().map(|&v| v as f64).collect();
        let x_nodes = DMatrix::from_row_slice(m, k, &data_vec);
        let x_final = build_mat_with_edges(&x_nodes, &self.selected_edges);

        let mut w_full_vec = Vec::new();
        if let Some(nw) = &self.node_weights {
            w_full_vec.extend(nw.iter());
        }
        if let Some(ew) = &self.edge_weights {
            w_full_vec.extend(ew.iter());
        }
        let w_full = DVector::from_vec(w_full_vec);

        let preds = (x_final * w_full).add_scalar(self.bias);

        Ok(PyArray1::from_vec(py, preds.data.as_vec().clone()))
    }

    fn get_node_weights<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        match &self.node_weights {
            Some(w) => Ok(PyArray1::from_vec(py, w.data.as_vec().clone())),
            None => Ok(PyArray1::from_vec(py, vec![])),
        }
    }

    fn get_edge_weights<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        match &self.edge_weights {
            Some(w) => Ok(PyArray1::from_vec(py, w.data.as_vec().clone())),
            None => Ok(PyArray1::from_vec(py, vec![])),
        }
    }

    fn get_selected_edges(&self) -> Vec<(usize, usize)> {
        self.selected_edges.clone()
    }

    fn dump(&self, path: &str) -> PyResult<()> {
        let serialized = bincode::serialize(self)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let mut file = File::create(path)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        file.write_all(&serialized)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    #[staticmethod]
    fn load(path: &str) -> PyResult<Self> {
        let mut file = File::open(path)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let mut serialized = Vec::new();
        file.read_to_end(&mut serialized)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let model = bincode::deserialize(&serialized)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(model)
    }
}

#[pymodule]
fn fast_ridge(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GraphModel>()?;
    Ok(())
}
