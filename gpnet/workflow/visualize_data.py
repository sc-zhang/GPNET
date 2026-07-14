from gpnet.io.data_io import GraphLoader, GraphSaver
from gpnet.io.message import Message
from gpnet.visualizer.visualizer import GraphVisualizer


def main(args):
    out_pic = args.output
    in_data = args.input
    node_cmap = args.node_cmap
    edge_cmap = args.edge_cmap
    node_size_ratio = args.node_ratio
    edge_width_ratio = args.edge_ratio
    is_predicted_only = args.predicted_only
    is_lower_better = args.lower
    highlight_color = args.highlight_color
    label_top_nodes = args.label_top_nodes
    label_top_edges = args.label_top_edges
    edge_alpha_min = args.edge_alpha_min
    draw_top_nodes = args.draw_top_nodes
    draw_top_edges = args.draw_top_edges

    Message.info("Loading data")
    graph_data = GraphLoader()
    graph_data.load_data(
        in_data,
        is_lower_better,
        is_predicted_only,
    )

    Message.info("Plotting")
    gv = GraphVisualizer(
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
    )
    gv.draw()

    Message.info("Saving data")
    out_pre = ".".join(out_pic.split(".")[:-1])
    ds = GraphSaver(out_pre)
    ds.save_graph(graph_data)
    Message.info("Finished")
