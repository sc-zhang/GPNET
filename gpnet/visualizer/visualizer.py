import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from gpnet.io.message import Message

mpl.use("Agg")


class GraphVisualizer(object):
    def __init__(
        self,
        out_pic,
        graph_data,
        node_size_ratio,
        edge_width_ratio,
        is_lower_better,
        highlight_color,
        node_cmap,
        edge_cmap,
        label_top_nodes,
        label_top_edges,
        edge_alpha_min,
        draw_top_nodes,
        draw_top_edges,
    ):
        self.__out_pic = out_pic
        self.__graph_data = graph_data
        self.__node_size_ratio = node_size_ratio
        self.__edge_width_ratio = edge_width_ratio
        self.__is_lower_better = is_lower_better
        self.__highlight_color = highlight_color
        self.__node_cmap = node_cmap
        self.__edge_cmap = edge_cmap
        self.__label_top_nodes = label_top_nodes
        self.__label_top_edges = label_top_edges
        self.__edge_alpha_min = edge_alpha_min
        self.__draw_top_nodes = draw_top_nodes
        self.__draw_top_edges = draw_top_edges

    @staticmethod
    def __gen_mapper(vmin, vmax, cmap_name, cmap_parts=100):
        if vmin == vmax:
            vmin = 0
            vmax = 100
        cmap = plt.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        mapper.set_array(np.arange(vmin, vmax, (vmax - vmin) * 1.0 / cmap_parts))
        return mapper

    @staticmethod
    def __resolve_node_overlaps(pos, graph, node_sizes, axes, iterations=50):
        figure = axes.figure
        bbox = axes.get_position()
        fig_width_in, fig_height_in = figure.get_size_inches()

        ax_width_in = bbox.width * fig_width_in
        ax_height_in = bbox.height * fig_height_in

        x_lim = axes.get_xlim()
        y_lim = axes.get_ylim()
        data_width = x_lim[1] - x_lim[0]
        data_height = y_lim[1] - y_lim[0]

        pts_per_data_x = (ax_width_in * 72.0) / data_width if data_width > 0 else 1.0
        pts_per_data_y = (ax_height_in * 72.0) / data_height if data_height > 0 else 1.0

        node_radius_data = {}
        for node in graph.nodes():
            radius_pt = np.sqrt(node_sizes[node] / np.pi)
            rx = radius_pt / pts_per_data_x
            ry = radius_pt / pts_per_data_y
            node_radius_data[node] = (rx + ry) / 2.0

        node_list = list(graph.nodes())
        for _ in range(iterations):
            any_moved = False
            for i in range(len(node_list)):
                for j in range(i + 1, len(node_list)):
                    ni, nj = node_list[i], node_list[j]
                    dx = pos[ni][0] - pos[nj][0]
                    dy = pos[ni][1] - pos[nj][1]
                    dist = np.sqrt(dx * dx + dy * dy)
                    if dist < 1e-9:
                        dx = (np.random.random() - 0.5) * 0.01
                        dy = (np.random.random() - 0.5) * 0.01
                        dist = np.sqrt(dx * dx + dy * dy)
                    min_dist = node_radius_data[ni] + node_radius_data[nj]
                    if dist < min_dist:
                        overlap = (min_dist - dist) / 2.0
                        px = (dx / dist) * overlap
                        py = (dy / dist) * overlap
                        pos[ni] = (pos[ni][0] + px, pos[ni][1] + py)
                        pos[nj] = (pos[nj][0] - px, pos[nj][1] - py)
                        any_moved = True
            if not any_moved:
                break

    def draw(self):
        has_label_n = self.__label_top_nodes > 0
        has_label_edges = self.__label_top_edges > 0

        if not has_label_n and has_label_edges:
            self.__label_top_nodes = self.__label_top_edges
            has_label_n = True
        if not has_label_edges and has_label_n:
            self.__label_top_edges = self.__label_top_nodes
            has_label_edges = True

        show_labels = has_label_n or has_label_edges

        node_weights, edge_weights = (
            self.__graph_data.node_weights,
            self.__graph_data.edge_weights,
        )

        valid_edges = [
            (g1, g2, w)
            for g1, g2, w in edge_weights
            if g1 in node_weights and g2 in node_weights
        ]

        rank_reverse = not self.__is_lower_better

        sorted_nodes = sorted(
            node_weights.keys(), key=lambda n: node_weights[n], reverse=rank_reverse
        )
        node_number_map = {}
        for idx, node in enumerate(sorted_nodes, start=1):
            node_number_map[node] = idx

        edges_sorted = sorted(valid_edges, key=lambda e: e[2])
        edges_sorted_by_weight = sorted(
            valid_edges, key=lambda e: e[2], reverse=rank_reverse
        )

        top_node_list = []
        top_edge_list = []

        if has_label_n:
            top_node_list = sorted_nodes[: self.__label_top_nodes]

        if has_label_edges and edges_sorted_by_weight:
            top_edge_list = [
                (u, v) for u, v, _ in edges_sorted_by_weight[: self.__label_top_edges]
            ]

        if self.__draw_top_edges < 0 and self.__draw_top_nodes > 0:
            self.__draw_top_edges = min(self.__draw_top_nodes, len(valid_edges))

        filter_active = (
            self.__draw_top_nodes > 0 and self.__draw_top_nodes < len(node_weights)
        ) or (self.__draw_top_edges > 0 and self.__draw_top_edges < len(valid_edges))

        if filter_active:
            priority_nodes = set(top_node_list)
            for u, v in top_edge_list:
                priority_nodes.add(u)
                priority_nodes.add(v)

            draw_edge_count = (
                min(self.__draw_top_edges, len(edges_sorted_by_weight))
                if self.__draw_top_edges > 0
                else len(edges_sorted_by_weight)
            )
            selected_edges = edges_sorted_by_weight[:draw_edge_count]

            for u, v, _ in selected_edges:
                priority_nodes.add(u)
                priority_nodes.add(v)

            drawn_nodes = set(priority_nodes)

            if self.__draw_top_nodes > 0 and len(drawn_nodes) < self.__draw_top_nodes:
                for node in sorted_nodes:
                    if node in drawn_nodes:
                        continue
                    if len(drawn_nodes) >= self.__draw_top_nodes:
                        break
                    drawn_nodes.add(node)

            final_edges = [
                (u, v, w)
                for u, v, w in selected_edges
                if u in drawn_nodes and v in drawn_nodes
            ]
        else:
            drawn_nodes = set(node_weights.keys())
            final_edges = edges_sorted_by_weight

        top_node_list = [n for n in top_node_list if n in drawn_nodes]
        top_edge_list = [
            (u, v) for u, v in top_edge_list if u in drawn_nodes and v in drawn_nodes
        ]

        labeled_nodes = set()
        node_label_map = {}
        for node in top_node_list:
            labeled_nodes.add(node)
        for u, v in top_edge_list:
            labeled_nodes.add(u)
            labeled_nodes.add(v)
        for node in labeled_nodes:
            node_label_map[node] = str(node_number_map[node])

        graph = nx.Graph()
        for node in drawn_nodes:
            graph.add_node(node, weight=node_weights[node])
        for u, v, w in final_edges:
            graph.add_edge(u, v, weight=w)

        num_nodes = graph.number_of_nodes()
        num_edges = graph.number_of_edges()

        Message.info("\t%d nodes, %d edges added" % (num_nodes, num_edges))

        node_weights_raw = [graph.nodes[n]["weight"] for n in graph.nodes()]
        if self.__is_lower_better:
            node_weights_raw = [-w for w in node_weights_raw]
        node_weight_min = min(node_weights_raw) if node_weights_raw else 0
        node_weight_max = max(node_weights_raw) if node_weights_raw else 1
        node_weight_range = node_weight_max - node_weight_min
        if node_weight_range == 0:
            node_weight_range = 1
        node_size_min_frac = 0.04

        node_sizes = {
            node: self.__node_size_ratio
            * np.sqrt(
                node_size_min_frac
                + (1.0 - node_size_min_frac)
                * ((node_weights_raw[i] - node_weight_min) / node_weight_range)
            )
            for i, node in enumerate(graph.nodes())
        }

        edge_abs_weights = [abs(d["weight"]) for _, _, d in graph.edges(data=True)]
        max_edge_abs = max(edge_abs_weights) if edge_abs_weights else 1

        edge_widths_by_edge = {}
        edge_alpha_by_edge = {}
        alpha_ranks = {}
        if edge_abs_weights:
            sorted_graph_edges = sorted(
                graph.edges(data=True), key=lambda e: e[2]["weight"]
            )
            for rank, (u, v, d) in enumerate(sorted_graph_edges):
                edge_widths_by_edge[(u, v)] = (
                    self.__edge_width_ratio * abs(d["weight"]) / max_edge_abs
                )
                alpha_rank = rank / max(len(sorted_graph_edges) - 1, 1)
                edge_alpha_by_edge[(u, v)] = (
                    self.__edge_alpha_min + (1.0 - self.__edge_alpha_min) * alpha_rank
                )
                alpha_ranks[(u, v)] = rank

        figure_base = min(20, max(9, num_nodes * 0.09))
        network_plot_size = figure_base

        layout_scale = 5.0 + 2.0 * np.sqrt(num_nodes)
        k_default = 6.0 / np.sqrt(num_nodes) if num_nodes > 0 else 1.0
        try:
            pos = nx.kamada_kawai_layout(graph)
        except Exception:
            pos = nx.spring_layout(
                graph, k=k_default, iterations=500, seed=42, scale=layout_scale
            )
        else:
            pos_array = np.array(list(pos.values()))
            center = pos_array.mean(axis=0)
            pos_array = (pos_array - center) * layout_scale
            for i, node in enumerate(pos):
                pos[node] = tuple(pos_array[i])

        figure = plt.figure(
            figsize=(
                network_plot_size * (1.5 if show_labels else 1.0),
                network_plot_size,
            ),
            dpi=300,
        )

        legend_ax = None
        if show_labels:
            gs = figure.add_gridspec(1, 2, width_ratios=[6, 1], wspace=0.01)
            axes = figure.add_subplot(gs[0, 0])
            legend_ax = figure.add_subplot(gs[0, 1])
        else:
            axes = figure.add_subplot(1, 1, 1)

        axes.set_xlim(-layout_scale * 1.5, layout_scale * 1.5)
        axes.set_ylim(-layout_scale * 1.5, layout_scale * 1.5)
        axes.set_aspect("equal")

        self.__resolve_node_overlaps(pos, graph, node_sizes, axes, iterations=80)

        all_x = [p[0] for p in pos.values()]
        all_y = [p[1] for p in pos.values()]
        max_radius_data = 0.0
        if node_sizes:
            bbox = axes.get_position()
            fig_w_in, fig_h_in = figure.get_size_inches()
            ax_w_in = bbox.width * fig_w_in
            ax_h_in = bbox.height * fig_h_in
            x_lim_cur = axes.get_xlim()
            y_lim_cur = axes.get_ylim()
            data_w = x_lim_cur[1] - x_lim_cur[0]
            data_h = y_lim_cur[1] - y_lim_cur[0]
            pts_per_x = (ax_w_in * 72.0) / data_w if data_w > 0 else 1.0
            pts_per_y = (ax_h_in * 72.0) / data_h if data_h > 0 else 1.0
            for n in graph.nodes():
                r_pt = np.sqrt(node_sizes[n] / np.pi)
                rx_data = r_pt / pts_per_x
                ry_data = r_pt / pts_per_y
                radius_data = max(rx_data, ry_data)
                if radius_data > max_radius_data:
                    max_radius_data = radius_data

        x_min = min(all_x) - max_radius_data * 4.0 if all_x else -layout_scale
        x_max = max(all_x) + max_radius_data * 4.0 if all_x else layout_scale
        y_min = min(all_y) - max_radius_data * 4.0 if all_y else -layout_scale
        y_max = max(all_y) + max_radius_data * 4.0 if all_y else layout_scale

        x_range = x_max - x_min
        y_range = y_max - y_min
        max_range = max(x_range, y_range)
        x_center = (x_min + x_max) / 2.0
        y_center = (y_min + y_max) / 2.0

        margin = max_range * 0.06
        half = max_range / 2.0 + margin
        axes.set_xlim(x_center - half, x_center + half)
        axes.set_ylim(y_center - half, y_center + half)

        node_colors = [graph.nodes[n]["weight"] for n in graph.nodes()]
        vmax = max(abs(v) for v in node_colors) if node_colors else 1

        nodes_collection = nx.draw_networkx_nodes(
            graph,
            pos,
            node_size=[node_sizes[n] for n in graph.nodes()],
            node_color=node_colors,
            cmap=plt.get_cmap(self.__node_cmap),
            vmin=-vmax,
            vmax=vmax,
            alpha=0.95,
            edgecolors="#222222",
            linewidths=0.6,
            ax=axes,
        )
        nodes_collection.set_zorder(4)
        node_size_list = [node_sizes[n] for n in graph.nodes()]

        if num_edges > 0:
            edges_sorted = sorted(
                graph.edges(data=True), key=lambda e: alpha_ranks.get((e[0], e[1]), 0)
            )

            for u, v, d in edges_sorted:
                width = edge_widths_by_edge.get((u, v), 0)
                alpha = edge_alpha_by_edge.get((u, v), self.__edge_alpha_min)
                edge_color = d["weight"]
                edges_collection = nx.draw_networkx_edges(
                    graph,
                    pos,
                    edgelist=[(u, v)],
                    width=width,
                    edge_color=[edge_color],
                    edge_cmap=plt.get_cmap(self.__edge_cmap),
                    edge_vmin=-max_edge_abs,
                    edge_vmax=max_edge_abs,
                    alpha=alpha,
                    node_size=node_size_list,
                    arrows=True,
                    connectionstyle="arc3,rad=0.35",
                    ax=axes,
                )

                for ec in edges_collection:
                    ec.set_zorder(2)

        if show_labels:
            if top_node_list:
                highlight_nodes = nx.draw_networkx_nodes(
                    graph,
                    pos,
                    nodelist=top_node_list,
                    node_size=[node_sizes[n] for n in top_node_list],
                    node_color="none",
                    edgecolors=self.__highlight_color,
                    linewidths=2.2,
                    ax=axes,
                )
                highlight_nodes.set_zorder(5)

            if top_edge_list:
                for u, v in top_edge_list:
                    base_width = edge_widths_by_edge.get((u, v), 0)
                    highlight_edges = nx.draw_networkx_edges(
                        graph,
                        pos,
                        edgelist=[(u, v)],
                        width=base_width + 1.5,
                        edge_color=self.__highlight_color,
                        alpha=0.85,
                        node_size=node_size_list,
                        arrows=True,
                        connectionstyle="arc3,rad=0.35",
                        ax=axes,
                    )
                    for he in highlight_edges:
                        he.set_zorder(4)

        if node_label_map:
            base_font_size = max(8, min(14, 220 / np.sqrt(num_nodes)))
            label_pos_shifted = {k: (v[0], v[1] + 0.015) for k, v in pos.items()}
            norm = plt.Normalize(vmin=-vmax, vmax=vmax)
            for node, label_text in node_label_map.items():
                digits = len(label_text)
                font_size = base_font_size / (1.0 + 0.18 * (digits - 1))
                x, y = label_pos_shifted[node]
                fill_rgba = plt.get_cmap(self.__node_cmap)(
                    norm(graph.nodes[node]["weight"])
                )
                r, g, b = fill_rgba[0], fill_rgba[1], fill_rgba[2]
                luminance = 0.299 * r + 0.587 * g + 0.114 * b
                text_color = "white" if luminance < 0.55 else "black"
                axes.text(
                    x,
                    y,
                    label_text,
                    fontsize=font_size,
                    fontfamily="sans-serif",
                    fontweight="bold",
                    ha="center",
                    va="center",
                    color=text_color,
                    zorder=6,
                )

        node_mapper = self.__gen_mapper(-vmax, vmax, self.__node_cmap)
        node_cb = figure.colorbar(
            node_mapper,
            ax=axes,
            orientation="horizontal",
            fraction=0.035,
            shrink=0.25,
            pad=0.04,
        )
        node_cb.set_label("Node weight", fontsize=9, fontweight="normal")
        node_cb.ax.tick_params(labelsize=7)
        node_cb.outline.set_linewidth(0.3)

        edge_mapper = self.__gen_mapper(-max_edge_abs, max_edge_abs, self.__edge_cmap)
        edge_cb = figure.colorbar(
            edge_mapper,
            ax=axes,
            orientation="horizontal",
            fraction=0.035,
            shrink=0.25,
            pad=0.04,
        )
        edge_cb.set_label("Edge weight", fontsize=9, fontweight="normal")
        edge_cb.ax.tick_params(labelsize=7)
        edge_cb.outline.set_linewidth(0.3)

        axes.spines["top"].set_visible(False)
        axes.spines["right"].set_visible(False)
        axes.spines["left"].set_visible(False)
        axes.spines["bottom"].set_visible(False)
        axes.set_xticks([])
        axes.set_yticks([])

        axes.set_title(
            f"Gene Network\n" f"({num_nodes} nodes, {num_edges} edges)",
            fontsize=15,
            fontweight="normal",
            pad=14,
        )

        if show_labels and legend_ax is not None:
            legend_font_size = max(10, min(10, 190 / np.sqrt(len(labeled_nodes))))

            max_num_width = 0
            for node in labeled_nodes:
                w = len(str(node_number_map[node]))
                if w > max_num_width:
                    max_num_width = w
            num_fmt = f"{{:>{max_num_width}}}"

            parts = []
            if top_node_list:
                parts.append(f"Top-{len(top_node_list)} Nodes")
                for node in top_node_list:
                    num = node_number_map[node]
                    parts.append(f"{num_fmt.format(num)}: {node}")

            if top_edge_list:
                parts.append("")
                parts.append(f"Top-{len(top_edge_list)} Edges")
                for u, v in top_edge_list:
                    nu = node_number_map.get(u, "?")
                    nv = node_number_map.get(v, "?")
                    parts.append(
                        f"{num_fmt.format(nu)}<->{num_fmt.format(nv)}: {u}<->{v}"
                    )

            legend_text = "\n".join(parts)

            legend_ax.axis("off")
            legend_ax.text(
                0.02,
                0.5,
                legend_text,
                fontsize=legend_font_size,
                fontfamily="monospace",
                ha="left",
                va="center",
                transform=legend_ax.transAxes,
                linespacing=1.15,
            )

        figure.savefig(
            self.__out_pic,
            bbox_inches="tight",
            facecolor="white",
            pad_inches=0.3,
        )
        Message.info("\tPlot saved to: %s" % self.__out_pic)
        plt.close(figure)
