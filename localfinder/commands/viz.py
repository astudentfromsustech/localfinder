from localfinder.utils import get_plotly_default_colors, visualize_tracks


def main(args):
    # Extract arguments
    input_files = args.input_files
    output_file = args.output_file
    method = args.method
    region = args.region
    colors = args.colors

    # If colors are not provided, use default colors
    if not colors:
        colors = get_plotly_default_colors(len(input_files))

    # Parse region if provided
    region_tuple = None
    if region:
        chrom, start, end = region
        region_tuple = (chrom, int(start), int(end))

    visualize_tracks(
        input_files=input_files,
        output_file=output_file,
        method=method,
        region=region_tuple,
        colors=colors,
    )
