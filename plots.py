import pandas as pd
from plotnine import *

# ------------------------------------------------------ #
# Define  attributes for consistent styling across plots #
# ------------------------------------------------------ #
AXIS_LINE_COLOR = "#888686"
HV_LINE_COLOR   = "#e00202"  # A distinct red color for horizontal lines
HV_LINE_SIZE    = 1
HV_LINE_ALPHA   = 0.75


def plot_tsoppy_barplot(data         : pd.DataFrame = None,
                        x_var        : str          = None,
                        y_var        : str          = None,
                        fill_var     : str          = None,
                        guide_title  : str          = None,
                        x_lab        : str          = None,
                        y_lab        : str          = None,
                        cart_ylim    : tuple        = None,
                        title        : str          = None,
                        x_lab_angle  : int          = 45,
                        fig_size     : tuple        = (15, 6),
                        alpha_value  : float        = 0.8,
                        hline_y      : float        = None,
                        hline_alpha  : float        = HV_LINE_ALPHA,
                        hline_color  : str          = HV_LINE_COLOR,
                        hline_label  : str          = "",
                        hline_size   : float        = HV_LINE_SIZE,
                        ann_y_offset : float        = 0.0,
                        y_tick_step  : int          = None
                        ) -> ggplot:
    """
    Function to create a bar plot using plotnine (ggplot2 in Python) for TSOPPI run metrics.
    This function is designed to be flexible and can be used for various run metrics by passing the appropriate parameters.

    When hline_y is provided, a horizontal line will be added to the plot at the specified y-value, along with an annotation label if hline_label is provided.
    

    Parameters:
    - data: A pandas DataFrame containing the data to be plotted.
    - x_var: The name of the column in the DataFrame to be used for the x-axis.
    - y_var: The name of the column in the DataFrame to be used for the y-axis.
    - fill_var: The name of the column in the DataFrame to be used for filling the bars (e.g., different runs).
    - guide_title: The title for the legend guide corresponding to the fill variable.
    - x_lab: The label for the x-axis.
    - y_lab: The label for the y-axis.
    - cart_ylim: A tuple specifying the limits for the y-axis (e.g., (0, 104)).
    - title: The title of the plot.
    - x_lab_angle: The angle for the x-axis labels (default is 45 degrees).
    - fig_size: A tuple specifying the size of the figure (default is (15, 6)).
    - alpha_value: The transparency level for the bars (default is 0.8).
    - hline_y: The y-value at which to draw a horizontal line (optional).
    - hline_alpha: The transparency level for the horizontal line (default is 1.0).
    - hline_color: The color of the horizontal line (default is "red").
    - hline_label: The label for the horizontal line annotation (default is an empty string).
    - hline_size: The thickness of the horizontal line (default is 1.0).
    - ann_y_offset: The vertical offset for the annotation label (default is 0.0, which means it will be placed directly at the end of the line). Negative values will place the label below the line, while positive values will place it above the line.
    - y_tick_step: Optional fixed separation between y-axis ticks (e.g., 10, 50, 100).

    Returns:
    - A plotnine ggplot object representing the bar plot. 
    """

    plot = (
        ggplot(data, aes(x = x_var, y = y_var)) +
        geom_col(aes(fill = fill_var), alpha = alpha_value) +
        #scale_fill_brewer(type = 'qual', palette = 'Set3') +
        guides(fill = guide_legend(guide_title)) +
        xlab(x_lab) +
        ylab(y_lab) +
        coord_cartesian(ylim = cart_ylim) +
        ggtitle(title) +
        theme_minimal() +
        theme(axis_text_x        = element_text(angle = x_lab_angle, size = 7, ha = 'right', hjust = 1.5),
              axis_ticks_major_x = element_line(color = AXIS_LINE_COLOR),
              figure_size        = fig_size,
              axis_line          = element_line(color = AXIS_LINE_COLOR),
              plot_title         = element_text(weight = 'bold', ha = 'left')
            ) +
        scale_x_discrete(limits = data[x_var].drop_duplicates().tolist(),
                         breaks = data[x_var].drop_duplicates().tolist(),
                         expand = (0.05, 0, 0.1, 0))
    )

    # Add fixed y-axis tick spacing when requested
    if y_tick_step is not None and y_tick_step > 0:
        y_series = pd.to_numeric(data[y_var], errors='coerce')
        data_max = y_series.max()
        upper_lim = cart_ylim[1] if cart_ylim is not None and cart_ylim[1] is not None else data_max

        if pd.notna(upper_lim):
            upper_tick = int(upper_lim) + y_tick_step
            y_breaks = list(range(0, upper_tick + 1, y_tick_step))
            plot = plot + scale_y_continuous(breaks = y_breaks)

    # Add horizontal line and annotation if hline_y is provided
    # Also add extra space on the right for the annotation label
    if hline_y is not None:
        plot = (plot +
                geom_hline(yintercept = hline_y,
                           alpha      = hline_alpha,
                           color      = hline_color,
                           size       = hline_size) +
                annotate("text", len(data) + 1.25, hline_y + ann_y_offset,
                         label = hline_label, angle = 90, size = 12)
               )   
    
    return plot



def plot_contamination_scatter(data                    : pd.DataFrame = None,
                               color_var               : str          = None,
                               label_var               : str          = None,
                               guide_title             : str          = 'Run',
                               title                   : str          = None,
                               x_lab_angle             : int          = 45,
                               fig_size                : tuple        = (15, 6),
                               max_contamination_score : float        = 5000,
                               usl_contamination_score : float        = None,
                               usl_contamination_pval  : float        = None
                               ) -> ggplot:
    """
    Function to create a contamination assessment scatter plot for TSOPPI run metrics.
    
    Parameters:
    - data: A pandas DataFrame containing the contamination data to be plotted.
    - color_var: The name of the column to use for point colors (e.g., 'RUN' or 'highlighted_run').
    - label_var: The name of the column to use for point labels (optional, set to None for no labels).
    - guide_title: The title for the legend guide.
    - title: The title of the plot.
    - x_lab_angle: The angle for the x-axis labels (default is 45 degrees).
    - fig_size: A tuple specifying the size of the figure (default is (15, 6)).
    - max_contamination_score: The maximum contamination score for x-axis range.
    - usl_contamination_score: The USL guideline value for contamination score (x-axis).
    - usl_contamination_pval: The USL guideline value for contamination p-value (y-axis).
    
    Returns:
    - A plotnine ggplot object representing the contamination scatter plot.
    """
    
    # Calculate breaks for axes
    x_breaks = [i * 1000 for i in range(int(max_contamination_score / 1000) + 2)]
    y_breaks = [i * 0.1 for i in range(11)]
    
    plot = (
        ggplot(data, aes(x='DNA_CONTAMINATION_SCORE', y='DNA_CONTAMINATION_P_VALUE')) +
        geom_point(aes(color=color_var), size=3) +
        guides(color=guide_legend(guide_title)) +
        xlab('Contamination score') +
        ylab('Contamination P-value') +
        coord_cartesian(xlim=(0, max_contamination_score + 500), ylim=(0.0, 1.05)) +
        scale_x_continuous(breaks=x_breaks, minor_breaks=[i * 500 for i in range(int(max_contamination_score / 500) + 2)]) +
        scale_y_continuous(breaks=y_breaks, minor_breaks=[i * 0.05 for i in range(21)]) +
        ggtitle(title) +
        annotate('rect', 
                 xmin=usl_contamination_score, 
                 xmax=max_contamination_score + 250,
                 ymin=usl_contamination_pval, 
                 ymax=1, 
                 fill=HV_LINE_COLOR, 
                 alpha=0.3) +
        geom_hline(yintercept=usl_contamination_pval, alpha=HV_LINE_ALPHA, color=HV_LINE_COLOR, size = HV_LINE_SIZE) +
        annotate('text', 
                 x=max_contamination_score + 500, 
                 y=usl_contamination_pval + 0.1,
                 label=f'USL_Guideline: {usl_contamination_pval}', 
                 angle=90, 
                 size=10) +
        geom_vline(xintercept=usl_contamination_score, alpha=HV_LINE_ALPHA, color=HV_LINE_COLOR, size = HV_LINE_SIZE) +
        annotate('text', 
                 x=usl_contamination_score + 200, 
                 y=1.05,
                 label=f'USL_Guideline: {usl_contamination_score}', 
                 size=10) +
        theme_minimal() +
        theme(
            axis_text_x=element_text(angle=x_lab_angle, vjust=0.5, hjust=1, size=7),
            figure_size=fig_size,
            axis_line=element_line(color=AXIS_LINE_COLOR),
            plot_title=element_text(weight='bold', ha='left'),
            panel_grid_major=element_line(color=AXIS_LINE_COLOR, size=0.75, alpha=0.5),
            panel_grid_minor=element_line(color=AXIS_LINE_COLOR, size=0.25, alpha=0.5)
        )
    )
    
    # Add text labels if label_var is provided
    if label_var is not None:
        plot = plot + geom_text(aes(label=label_var), angle=-25, nudge_y=0.02, size=7)
    
    return plot