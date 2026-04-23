import pandas as pd
from plotnine import *


# ggplot(plot_data3, aes(x='RUN_ID', y='PCT_PF_READS')) +
# geom_col(aes(fill='RUN_ID'), alpha=0.8) +
# guides(fill=guide_legend('Run')) +
# xlab('\\nRun ID') +
# ylab('Percentage of PF reads') +
# coord_cartesian(ylim=(0, 104)) +
# ggtitle('[LocalApp core run metric] Percentage of reads passing instrument filters') +
# geom_hline(yintercept=80, alpha=0.3, color='red') +
# annotate('text', x=joint_qc_table.height + 1, y=80, label='LSL_Guideline: 80', angle=90) +
# theme_minimal() +
# theme(axis_text_x=element_text(angle=angle_x_names, vjust=0.5, hjust=1, size=7),
#       figure_size=(15, 6),
#       axis_line=element_line(color="#4d4d4d"),
#       plot_title=element_text(weight='bold'))


# ggplot(plot_data, aes(x='RUN_ID', y='CLUSTERS_PASSING_FILTER')) +
# geom_col(aes(fill='RUN_ID'), alpha=0.8) +
# guides(fill=guide_legend('Run')) +
# xlab('Run ID') +
# ylab('Clusters passing filter (percent)') +
# coord_cartesian(ylim=(0, 104)) +
# ggtitle('[Sequencer flowcell run metric] Clusters passing filter '
#         '(based on RunCompletionStatus.xml files)\n'
#         '(only displaying runs with non-NA values)') +
# theme_minimal() +
# theme(axis_text_x=element_text(angle=angle_x_names, vjust=0.5, hjust=1, size=7),
#       figure_size=(15, 6),
#       axis_line=element_line(color="#4d4d4d"),
#       plot_title=element_text(weight='bold'))


def plot_tsoppy_barplot(data        : pd.DataFrame = None,
                        x_var       : str          = None,
                        y_var       : str          = None,
                        fill_var    : str          = None,
                        guide_title : str          = None,
                        x_lab       : str          = None,
                        y_lab       : str          = None,
                        cart_ylim   : tuple        = None,
                        title       : str          = None,
                        x_lab_angle : int          = 45,
                        fig_size    : tuple        = (15, 6),
                        alpha_value : float        = 0.8,
                        hline_y     : float        = None,
                        hline_alpha : float        = 1.0,
                        hline_color : str          = "red",
                        hline_label : str          = "",
                        hline_size  : float        = 1.0
                        ) -> ggplot:
    """
    Function to create a bar plot using plotnine (ggplot2 in Python) for TSOPPI run metrics.
    This function is designed to be flexible and can be used for various run metrics by passing the appropriate parameters.
    
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
        theme(axis_text_x = element_text(angle = x_lab_angle, vjust = 0.5, hjust = 1, size = 7),
              figure_size = fig_size,
              axis_line   = element_line(color = "#4d4d4d"),
              plot_title  = element_text(weight = 'bold')
            )
    )

    # Add horizontal line and annotation if hline_y is provided    
    if hline_y is not None:
        plot = (plot +
                geom_hline(yintercept = hline_y,
                           alpha      = hline_alpha,
                           color      = hline_color,
                           size       = hline_size) +
                annotate("text", len(data) + 0.49, hline_y, label = hline_label, angle = 90, size = 13)
               )   
    
    return plot