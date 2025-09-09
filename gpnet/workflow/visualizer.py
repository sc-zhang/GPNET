from gpnet.io.data_io import GraphLoader, GraphSaver
from gpnet.io.message import Message
from pyecharts import options as opts
from pyecharts.charts import Graph


def plot_chart(nodes, edges, categories, out_html):
    g = (
        Graph(init_opts={"width": "1080px", "height": "1080px"})
        .add(
            "",
            nodes,
            edges,
            categories,
            repulsion=500,
            layout="circular",
            is_rotate_label=True,
            linestyle_opts=opts.LineStyleOpts(curve=0.2),
            label_opts=opts.LabelOpts(is_show=False)
        )
        .set_global_opts(
            legend_opts=opts.LegendOpts(is_show=False),
            title_opts=opts.TitleOpts(title="")
        )
        .render(out_html)
    )


def main(args):
    out_html = args.output
    in_data = args.input
    node_cmap = args.node_cmap
    edge_cmap = args.edge_cmap
    node_size_ratio = args.node_ratio
    edge_width_ratio = args.edge_ratio
    is_lower_better = args.lower

    Message.info("Loading data")
    dl = GraphLoader()
    dl.load_data(in_data, is_lower_better, node_size_ratio, edge_width_ratio, node_cmap, edge_cmap)

    Message.info("Plotting")
    plot_chart(dl.nodes, dl.edges, dl.categories, out_html)

    out_pre = out_html if not out_html.endswith(".html") else '.'.join(out_html.split('.')[:-1])
    ds = GraphSaver(out_pre)
    ds.save_graph(dl)
    Message.info("Finished")
