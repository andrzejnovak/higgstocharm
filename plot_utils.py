def format_legend(ax, ncols=2, handles_labels=None, **kwargs):
    if handles_labels is None:
        handles, labels = ax.get_legend_handles_labels()
    else:
        handles, labels = handles_labels
    nentries = len(handles)

    kw = dict(framealpha=1, **kwargs)
    split = nentries // ncols * ncols
    leg1 = ax.legend(handles=handles[:split],
                     labels=labels[:split],
                     ncol=ncols,
                     loc="upper right",
                     **kw)
    if nentries % 2 == 0:
        return leg1

    ax.add_artist(leg1)
    leg2 = ax.legend(handles=handles[split:],
                     labels=labels[split:],
                     ncol=nentries - nentries // ncols * ncols,
                     **kw)

    leg2.remove()

    leg1._legend_box._children.append(leg2._legend_handle_box)
    leg1._legend_box.stale = True
    return leg1
