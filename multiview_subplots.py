def multiview_subplots(fig_params, plot_params):
    """
    Create multiview whole-brain surface subplots with a common colorbar displayed
    
    params:
    fig_params (dict): Dictionary containing figure parameters
    
    plot_params (dict): Dictionary containing parameters relating to plotting actual data
    """
    
    plt.style.use('seaborn-poster')
    
    # Create global colorbar
    show_cmap = plot_params['show_cmap']
    cmap = plot_params['cmap']
    data = plot_params['data']
    threshold = plot_params['threshold']
    vmax = plot_params['vmax']
    cbar_label = plot_params['cbar_label']
    
    # Create subplots
    dims = fig_params['dim']
    figsize = fig_params['figsize']
    title = fig_params['title']
    fname = fig_params['fname']
    surf_mesh = plot_params['surf_mesh']
    bg_map = plot_params['bg_map']
    hemis = plot_params['hemis']
    views = plot_params['views']
    symmetric = plot_params['symmetric']
    
    fig, axes = plt.subplots(nrows=dims[0], ncols=dims[1], subplot_kw={'projection': '3d'}, figsize=figsize)
                             
    ax_flat = axes.flat

    vmin = threshold
    vmax = vmax
    
    our_cmap = get_cmap(cmap)
    norm = Normalize(vmin=vmin, vmax=vmax)

    for i, ax in enumerate(ax_flat):
        nilearn.plotting.plot_surf(
            surf_mesh=surf_mesh[i],
            surf_map=data[:int(len(data) / 2)] if hemis[i] == 'left' else data[int(len(data) / 2):],
            hemi=hemis[i],
            view=views[i], 
            bg_map=bg_map[i], cbar_vmin=vmin, cbar_vmax=vmax, cmap=cmap, 
            symmetric_cmap=symmetric, vmin=vmin, vmax=vmax, threshold=threshold, axes=ax
        )
                             
        ax.dist = 5.5  # Decrease to zoom-in
    if show_cmap:
        proxy_mappable = ScalarMappable(cmap=our_cmap, norm=norm)
        proxy_mappable.set_array(data)

        nb_ticks = 2
        ticks = np.linspace(vmin, vmax, nb_ticks)
        bounds = np.linspace(vmin, vmax, our_cmap.N)

        p0 = ax_flat[len(ax_flat) - (dims[1] - 1)].get_position().get_points().flatten()
        p1 = ax_flat[len(ax_flat) - 2].get_position().get_points().flatten()
        ax_cbar = fig.add_axes([p0[0], 0, p1[2]-p0[0], 0.025])
        plt.colorbar(proxy_mappable, cax=ax_cbar, boundaries=bounds, ticks=ticks, spacing='proportional', orientation='horizontal', pad=0, format='%.0e')
        ax_cbar.set_xlabel(cbar_label, fontsize=28)
    
    fig.text(0, 0.71, title, fontsize=28)
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0)
    if fname:
        plt.savefig(fname, bbox_inches='tight', pad_inches=0.05)
