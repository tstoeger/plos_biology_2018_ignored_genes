import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool

# __all__ = ['colored_bokeh']


def colored_bokeh(
        data,
        tooltips,
        output_path,
        minimum=None,
        maximum=None,
        colormap='inferno',
        output_name='',
        circle_size=3,
        alpha=0.5,
):
    """
    #data=dict(
    #            x = Y[:, 0],
    #            y = Y[:, 1],
    #            z = Z,
    #            desc=glabel,
    #        )

    #tooltips=[
    #            ("desc", "@desc"),
    #            ("att", "@knowledge"),
    #        ]
    """

    Z = data['z']       # will be displayed as color

    if minimum is None:
        minimum = min(Z)

    if maximum is None:
        maximum = max(Z)

    norm = mpl.colors.Normalize(vmin=minimum, vmax=maximum)

    cmap = plt.get_cmap(colormap)

    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    colors = []
    for j in Z:
        v = m.to_rgba(j)
        r = "%0.2X" % int(round(v[0] * 255))
        g = "%0.2X" % int(round(v[1] * 255))
        b = "%0.2X" % int(round(v[2] * 255))
        s = '#' + r + g + b
        colors.append(s)

    # o = 'bokeh_' + output_name + '.html'
    output_file(output_path)

    source = ColumnDataSource(
        data
    )

    hover = HoverTool(
        tooltips=tooltips
    )

    p = figure(plot_width=1400, plot_height=1400, tools=[hover],
               title=output_name)

    p.circle('x', 'y', size=circle_size, fill_color=colors,
             source=source, line_color=None, fill_alpha=alpha)

    show(p)
