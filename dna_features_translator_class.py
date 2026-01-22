from dna_features_viewer import BiopythonTranslator
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import annotation
from annotation import record_check
from Bio import SeqIO
from os import getcwd
from sys import platform
import numpy as np 
import os
import utils
from mycolorpy import colorlist as mcp
from importlib import reload
from scale_axis import AsymScale
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib import patches
mscale.register_scale(AsymScale)

reload(annotation)
reload(utils)
PATH = getcwd()

if platform != 'win32':
    slash = '/'
else:
    slash = '\\'

plt.style.use('ggplot')

class MitoTranslator(BiopythonTranslator):
  """
  A class of custom coloration and filtering for usage by dna_features_viewer
  """
  def compute_feature_color(self, feature):
    """
    Colors CDS blue, misc_feature grey, rRNA green and introns red
    otherwise gold
    """
    if feature.type == "CDS":
      # Complex I : Blue
      # Complex III : Teal #03fcf8
      # Complex IV : Magenta #f803f8
      # Complex V : Pink #f005d0
      try:
        gene = utils.Which_tRNA(feature.qualifiers['product'], verbosity = 0)
        if 'ND' in gene.upper() or 'nad' in gene:
          return 'blue'
        elif 'COX' in gene.upper():
          return '#03fcf8'
        elif 'CYT' in gene.upper() or 'cob' in gene:
          return '#c1ffb0'
        elif 'ATP' in gene.upper():
          return '#f005d0'
        else:
          return 'orange'
      except (KeyError, ValueError): return 'orange'
    elif feature.type == 'misc_feature':
      return 'grey'
    elif feature.type == 'rRNA':
      return 'green'
    elif feature.type == 'intron':
      return 'red'
    elif feature.type == 'tRNA':
      return 'purple'
    elif feature.type == 'promoter':
      return 'black'
    else:
      return 'gold'

  def compute_feature_label(self, feature):
    """
    If feature is not CDS misc rRNA or intron, dont label it
    """
    wanted = ['rRNA','intron','gene', 'CDS', 'promoter']

    if feature.type not in wanted:
      return None
    else:
      feature.qualifiers['product'] = utils.Which_tRNA(feature.qualifiers['product'][0], verbosity = 0)
      try: feature.qualifiers['gene'] = utils.Which_tRNA(feature.qualifiers['gene'][0], verbosity = 0)
      except KeyError: pass
      return BiopythonTranslator.compute_feature_label(self, feature)

  def compute_filtered_features(self, features):
    """
    If feature not in wanted, dont graph it aswell"""
    wanted = ['CDS','misc_feature','rRNA','intron', 'tRNA']
    return [feature for feature in features if feature.type in wanted]

class MitoTranslatorT(BiopythonTranslator):
  """
  A class of custom coloration and filtering for usage by dna_features_viewer
  """
  def compute_feature_color(self, feature):
      """
      Colors CDS blue, misc_feature grey, rRNA green and introns red
      otherwise gold
      """
      if feature.type == "CDS":
        # Complex I : Blue
        # Complex III : Teal #03fcf8
        # Complex IV : Magenta #f803f8
        # Complex V : Pink #f005d0
        try:
          gene = utils.Which_tRNA(feature.qualifiers['product'][0], verbosity = 0)
          if 'ND' in gene.upper() or 'nad' in gene:
            return 'blue'
          elif 'COX' in gene.upper():
            return '#03fcf8'
          elif 'CYT' in gene.upper() or 'cob' in gene:
            return '#c1ffb0'
          elif 'ATP' in gene.upper():
            return '#f005d0'
          else:
            return 'orange'
        except (KeyError, ValueError): return 'orange'
      elif feature.type == 'misc_feature':
        return 'grey'
      elif feature.type == 'rRNA':
        return 'green'
      elif feature.type == 'intron':
        return 'red'
      elif feature.type == 'tRNA':
        return 'purple'
      elif feature.type == 'promoter':
        return 'black'
      else:
        return 'gold'

  def compute_feature_label(self, feature):
    """
    If feature is not CDS misc rRNA or intron, dont label it
    """
    wanted = ['rRNA','intron','gene', 'CDS', 'tRNA']
    if feature.type not in wanted:
      return None
    else:
      feature.qualifiers['product'] = utils.Which_tRNA(feature.qualifiers['product'][0], verbosity = 0)
      if feature.qualifiers['product'] == None:
        try: feature.qualifiers['gene'] = utils.Which_tRNA(feature.qualifiers['gene'][0], verbosity = 0)
        except KeyError: pass
      else:
        feature.qualifiers['gene'] = feature.qualifiers['product']
      return BiopythonTranslator.compute_feature_label(self, feature)

  def compute_filtered_features(self, features):
    """
    If feature not in wanted, dont graph it aswell"""
    wanted = ['CDS','misc_feature','rRNA','intron', 'tRNA']
    return [feature for feature in features if feature.type in wanted]
    
def modify_y(y, plot, log, to_window, window):
  """
  Helper function to coverage_graph. Takes the y-axis and modifies it according to the parameters sent to coverage_graph.

  Parameters
  ----------
  y : list
    List of y-axis values
  plot : str
    Name of the df column to be plotted
  log : bool
    If True, logarithmic y-axis
  
  Returns
  -------
  y : list
    List of y-axis values (modified)
  """
  if any(i in plot for i in to_window):
    y = y.rolling(window = window, min_periods = 1).mean()
  if log:
    y = np.log10(y)
  return y

def make_cols_neg(df):
  """
  Helper function to make all columns negative
  """
  df = df.copy()
  for i in df.columns:
    if 'neg' in i and df[i].max() > 0:
      df[i] = df[i] * -1
  return df

def plot_bool(df, col, ax, addit_samples, colors, color):
  """
  Helper function to plot_coverage. Plots a boolean column (True/False) as a rectangle graph.
  """
  colors = [color] + colors
  df = df.copy()
  try:
    if addit_samples:
      if not isinstance(addit_samples, list): raise ValueError('addit_samples must be a list')
  except ValueError:
    addit_samples = [addit_samples]
  if addit_samples:
    n = len(addit_samples) + 1
    for i, _ in enumerate(addit_samples):
      addit_samples[i] = addit_samples[i] = addit_samples[i].loc[addit_samples[i].Position.isin(df.Position)]
    addit_samples = [df] + addit_samples
  else:
    n = 1
    addit_samples = [df]
  for i, _ in enumerate(addit_samples):
    # Check if the possible values of addit_samples[i][col] are 0 and 1 or 0 and -1
    values = [int(i) for i in addit_samples[i][col].unique()]
    if (-1 in values) and (0 in values):
      neg_mode = True
    else:
      neg_mode = False
    bed = utils.per_position_to_bed(addit_samples[i], reg_col = col)
    for interval in bed.itertuples():
      # Define rectangle
      if not neg_mode:
        rect = patches.Rectangle((interval.start, -1 * (1/n - 1/(i + 1))), interval.end - interval.start, 1/n, linewidth = 1, facecolor = colors[i], edgecolor = colors[i], alpha = 0.8)
        # Add the patch to the Axes
        ax.add_patch(rect)
      else: 
        rect = patches.Rectangle((interval.start,  (0 - 1/n)*(i + 1)), interval.end - interval.start, 1/n, linewidth = 1, facecolor = colors[i], edgecolor = colors[i], alpha = 0.8)
        # Add the patch to the Axes
        ax.add_patch(rect)       
  if neg_mode:
    # Set y range from -1 to 1
    ax.set_ylim(-1, 1)
  # Remove the y spine
  ax.spines['left'].set_visible(False)
  # Remove the y ticks
  ax.yaxis.set_ticks([])
  # Remove the yticklabels
  ax.set_yticklabels('')
    # Set the y axis to be invisible
  ax.set_ylabel(col.replace('neg_', '').replace('pos_', '').capitalize())
  ax.tick_params(axis = 'y', which = 'both', left = False, labelleft = True)
  return ax


def coverage_graph(
  org, sample, strand, addit_samples = None, labels = [''], plotlist = ['coverage'], prediction = None, transcripts = None, transcript_overlap_mode = False, window = 100, scale = 'linear', log = False, peaks = None, peaks_alpha = 1, peak_names = None, peaks_annotate = False, peaks_color = 'black', tis_term_df = None, pos_color = 'red', neg_color = 'blue', neut_color = 'purple', to_window = ['z', 'coverage', 'RPM'], hline = '', title_org = True, custom_title = '', savefig = '', plot_range = None, alpha = 1, upscaling = 1, fill = False, y_up_lim = False, style = 'default', despine = True, return_ax = False, fig = None, axes = None, legend_outside = False, width = 15, ylabels = None, yfontsize = 12, tight_layout = False):
  """
  Create two graphs aligned, annotations and coverage.

  Parameters
  ----------
  org : str
    Organism name
  sample : pd.DataFrame
    sample df created by annotation.py's main
  strand : bool 
    strand of coverage graph
  addit_samples : pd.DataFrame, None or list
    additional sample df created by annotation.py's main
  labels : list
    list of labels for the samples
  plotlist : list
    A list of column names to plot on the y axis against position.
  prediction : pd.DataFrame
    A df of the prediction dataframe created by annotation.py's main
  transcripts : list
    nested list in this format [initiation, termination, strand] for the graph annotations.
  window : int
    rolling function value
  scale : str
    scale of the y-axis
  log : bool
    take the log values of the plotlist
  peaks : list
    A list of x axis values (positions) to insert a vertical line in
  peaks_alpha : float
    alpha value for the peaks
  peaks_annotate : bool
    If True, annotate the peaks
  peaks_color : str
    color of the peaks
  tis_term_df : pd.DataFrame
    A df of the tis_term dataframe created by annotation.py's main
  pos_color : str
    Positive expression coloration
  neg_color : str
    Negative expression coloration
  neut_color : str
    Neut expression coloration
  to_window : list
    List of sample df columns to apply the window function to
  hline : str
    If "mean" - include a horizontal dashed line at the mean value acrosss the y axis
  title_org : bool
    If True - include the organism name in the title
  custom_title : str
    If not empty - include a custom title 
  savefig : str
    If not empty - save the figure to this path
  plot_range : list
    If not empty - plot a range of values from the list
  alpha : float
    Transparency of the plot
  upscaling : float
    Upscaling factor for the plot
  fill : bool
    If True - fill the area under the plot
  
  Notes
  -----
  The plotlist can be a list of strings or a list of lists. If a list of lists, the first element of the list is the name of the plot and the second element is the name of the column in the sample df.
  Does not return anything, create graph inline
  """
  n = len(labels) if len(labels) > 5 else 5
  RANDOM_COLORS = mcp.gen_color('tab10', n = n)
  POS_COLORS = ['black', '#EADDCA', '#9300FF', '#00BA44', '#0093BA', '#0030BA', '#00abff', '#0000ff', '#ab00ff', '#ff00ff', '#ff00ab'] #colors for the positive strand
  NEG_COLORS = ['green', '#00e7ff', '#00ffd4', '#00ff9c', '#00ff6c', '#00ff3c', '#00ff00', '#3cff00', '#6cff00', '#9cff00', '#d4ff00'] #colors for the negative strand
  NEUT_COLORS = RANDOM_COLORS #colors for the neutral strand
  #TODO(Noam) - This 'kind of' works but is not very good. Need to go over line by line and fix this for multiple samples.
  plt.style.use(style)
  if strand == 'both':
    
    coverage_graph(org = org, sample = sample,strand = True, addit_samples = addit_samples, labels = labels, plotlist = plotlist, prediction = prediction,transcripts = transcripts,window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peak_names = peak_names,peaks_annotate =  peaks_annotate,peaks_color = peaks_color,  tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =  custom_title, savefig =  savefig[0: -4] + '_heavy' + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine, legend_outside = legend_outside, width = width, ylabels = ylabels, yfontsize = 12, tight_layout=tight_layout)
    coverage_graph(org = org, sample = sample,strand = False,addit_samples =  addit_samples,labels =  labels,  plotlist = plotlist, prediction = prediction,transcripts = transcripts,window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peak_names = peak_names,peaks_annotate =  peaks_annotate,peaks_color = peaks_color,  tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =   custom_title,savefig =  savefig[0: -4] + '_light' + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine, legend_outside = legend_outside, width = width, ylabels = ylabels, yfontsize = 12, tight_layout=tight_layout)
    return
  
  elif strand == 'overlap' or strand == 'overlap_noneg':
    fig, axes = coverage_graph(org = org, sample = sample,strand = True, addit_samples = addit_samples, labels = labels, plotlist = plotlist, prediction = prediction,transcripts = transcripts, transcript_overlap_mode = True, window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peak_names = peak_names,peaks_annotate =  peaks_annotate,peaks_color = peaks_color,  tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =  custom_title, savefig =  savefig[0: -4]  + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine, return_ax = True, legend_outside = legend_outside, width = width, ylabels = ylabels, yfontsize = 12, tight_layout=tight_layout)
    
    if isinstance(addit_samples, pd.DataFrame) or isinstance(addit_samples, list):
      if not isinstance(addit_samples, list): addit_samples = [addit_samples]
      if 'noneg' not in strand:
        for i, _ in enumerate(addit_samples):
          addit_samples[i] = make_cols_neg(addit_samples[i])
        sample = make_cols_neg(sample)

    coverage_graph(org = org, sample = sample, strand = False,addit_samples =  addit_samples,labels =  labels,  plotlist = plotlist, prediction = prediction,transcripts = transcripts, transcript_overlap_mode = True, window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peak_names = peak_names,peaks_annotate =  peaks_annotate, peaks_color = peaks_color, tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =   custom_title,savefig =  savefig[0: -4] + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine, axes = axes, fig = fig, return_ax = False, legend_outside = legend_outside, width = width, ylabels = ylabels, yfontsize = 12, tight_layout=tight_layout)
    axes[1].set_title('')
    axes[1].axhline(0, color = 'black', alpha = .8, linestyle = '-')
    for ax in axes[1:]:

      # Change the negative y axis ticks to positive labels only if ax.get_yticks() is not empty and not in the ranges 0 - 1 or -1 - 1
      ticks = [int(i) for i in ax.get_yticks()]
      if len(ticks) > 0 and not all([i in [0, 1] for i in ticks]) and not all([i in [-1, 0, 1] for i in ticks]):
        ax.set_yticks([i for i in ticks])
        ax.set_yticklabels([str(round(i)).replace('-','') for i in ticks])
    # Get handles and labels for the legend
    handles, labels = axes[1].get_legend_handles_labels()
    # Remove old legend
    axes[1].get_legend().remove()
    # Add heavy and light strand labels
    labels = [i + ' ' + j for i, j in zip(['Heavy strand', 'Heavy strand', 'Light strand', 'Light strand'], labels)]
    # Create a new legend
    if legend_outside:
      fig.legend(handles, labels, 
          ncol=2, fancybox=True, shadow=False)
    else:
      fig.legend(handles, labels, loc = 'upper right', bbox_to_anchor = (1, 1), bbox_transform = axes[1].transAxes)
    if savefig != None:
      plt.savefig(savefig, dpi = 300, bbox_inches = 'tight')

    return

  try: prediction = list(prediction)
  except TypeError: pass
  if plot_range == None:
    plot_range = [1, sample.Position.max()]
  sample = sample[sample.Position.isin(range(plot_range[0], plot_range[1] + 1))] # Limit the sample df to the plot range
  record_exists = record_check(org) # Check if the organism is in the saved and if not create it
  if not record_exists:
    gb_file_exists = record_check(org, mode = 'gb') # Check if the genbank file exists
    if not gb_file_exists:
      raise FileNotFoundError(f'{org} genbank file does not exist!')
    
  record = SeqIO.read(f'{PATH}{slash}genbank_DB{slash}{org}.gbk', format = 'genbank') # Read the genbank file
  if strand == False:
    plotlist = ['neg_' + i for i in plotlist]
    color = neg_color
    color_list = NEG_COLORS   
  elif strand == True:
    plotlist = ['pos_' + i for i in plotlist]
    color = pos_color
    color_list = POS_COLORS
  elif strand == 'both':
    plotlist_temp = ['pos_' + i for i in plotlist]
    plotlist = ['neg_' + i for i in plotlist]
    plotlist += plotlist_temp
  else:
    color = neut_color
    color_list = NEUT_COLORS
  if fig == None:
    fig, axes = plt.subplots(len(plotlist) + 1, 1, figsize = (width, upscaling * 1.5 * (len(plotlist) + 1)), sharex = True, gridspec_kw = {'height_ratios': [1/upscaling] + [1 for _ in range(len(plotlist))]})  # Create the figure and axes    
  x = sample.Position
  for i, plot in enumerate(plotlist):
    if despine:
      sns.despine(ax = axes[i + 1], top = True, right = True, trim = False)
    i+=1
    # If the sample[plot] series contains only 0 and 1, then transfer to plot_bool function and continue
    if len(sample[plot].unique()) == 2 and 0 in sample[plot].unique() and (1 in sample[plot].unique() or (-1 in sample[plot].unique())):
      plot_bool(sample, plot, axes[i], addit_samples, color_list, color)
      continue
    if strand == 'both':
      if 'pos_' in plot:
        color = pos_color
        y = sample[plot]
      elif 'neg_' in plot:
        color = neg_color
        y = sample[plot]
    else:
      try: y = sample[plot]
      except KeyError:
        try: 
          y = sample[plot.replace('neg_','').replace('pos_','')]
          print(f'Parameter {y} does not have a strand-specific representation!\n')
        except KeyError:
          print(f'Parameter {plot} does not exist!\n')
          continue
    y = modify_y(y, plot, log, to_window, window)
    sns.lineplot(ax = axes[i],x = x,y = y, color = color, label = labels[0], alpha = alpha, legend = False)
    if fill: axes[i].fill_between(x, y, color = color, alpha = alpha/2)
    try:
      if addit_samples:
        if not isinstance(addit_samples, list): raise ValueError('addit_samples must be a list')
    except ValueError:
      addit_samples = [addit_samples]
    if addit_samples:
      for j, addit_sample in enumerate(addit_samples):
        if len(addit_samples) + 1 > len(labels): raise ValueError('The number of labels must match the number of addit_samples')
        addit_sample = addit_sample[addit_sample.Position.isin(range(plot_range[0], plot_range[1] + 1))]
        try: y_cur = addit_sample[plot]
        except KeyError:
          try:
            y_cur = addit_sample[plot.replace('neg_','').replace('pos_','')]
          except KeyError:
            print(f'Skipping {plot} in sample {j + 1}')
            continue
        y_cur = modify_y(y_cur, plot, log, to_window, window)
        axes[i].plot(x, y_cur, color = color_list[j], label = labels[j+1], alpha = alpha)
        if fill: axes[i].fill_between(x, y_cur, color = color_list[j], alpha = alpha/2)

    if i == len(plotlist):
      axes[i].set_xlabel('Position')
      # Add more xticklabels
      plot_span = plot_range[1] - plot_range[0]
     # if plot_span < 1000:
     #   axes[i].set_xticks(np.arange(plot_range[0]-1, plot_range[1], 200))
     #   axes[i].set_xticklabels(np.arange(plot_range[0] - 1, plot_range[1], 200), rotation = 45)
      #else:
     #   axes[i].set_xticks(np.arange(plot_range[0] - 1, plot_range[1], 500))
     #   axes[i].set_xticklabels(np.arange(plot_range[0] - 1, plot_range[1], 500), rotation = 45)
    if not ylabels or len(ylabels) != len(plotlist):
      axes[i].set_ylabel(plot.replace('_', ' ').replace('neg ','').replace('pos ','').capitalize(), fontsize = yfontsize)
    else:
      axes[i].set_ylabel(ylabels[i - 1], fontsize = yfontsize)

    #If hline is requested, plot a horizontal line of the mean value
    if hline != '':
      if type(hline) == str:
        if hline == 'mean': # If the hline is the mean value
          if plot in to_window:
            y_pos = y[y > 1].mean()
          else:
            y_pos = y[y > 0].mean()
        elif hline == 'median': # If the hline is the median value
          if plot in to_window:
            y_pos = y[y > 1].median()
          else:
            y_pos = y[y > 0].median()
        elif hline == 'max': # If the hline is the max value
          if plot in to_window:
            y_pos = y[y > 1].max()
          else:
            y_pos = y[y > 0].max()
        elif hline == 'min': # If the hline is the min value
          if plot in to_window:
            y_pos = y[y > 1].min()
          else:
            y_pos = y[y > 0].min()
        elif '%' in hline: # If the hline is a percentage - return the requested quantile
          if plot in to_window:
            y_pos = y[y > 1].quantile(float(hline.replace('%', ''))/100)
          else:
            y_pos = y[y > 0].quantile(float(hline.replace('%', ''))/100)
      elif type(hline) == int:
        y_pos = hline
      else: y_pos = 1
      axes[i].axhline(y = y_pos, color = 'black', alpha = 0.7, linestyle = ':')
    
    #Plot terminations/initiations in the transcripts list which is in this format: [start,end,strand]
    if transcripts:
      c = 0
      # Normal mode
      if not transcript_overlap_mode:
        for t in transcripts:
          if t[2] == strand:
            if c == 0:
              axes[i].axvline(x = t[0], color = 'lawngreen', label = 'Initiation')
              axes[i].axvline(x = t[1], color = 'lightcoral', label = 'Termination')
              c+=1
            else:
              axes[i].axvline(x = t[0], color = 'lawngreen')
              axes[i].axvline(x = t[1], color = 'lightcoral')
          elif strand != True and strand != False:
            if c == 0:
              axes[i].axvline(x = t[0], color = 'lawngreen')
              axes[i].axvline(x = t[1], color = 'lightcoral')
              c+=1
            else:
              axes[i].axvline(x = t[0], color = 'lawngreen')
              axes[i].axvline(x = t[1], color = 'lightcoral')
      # Designed specifically for overlap mode
      else:
        c = 0
        if strand != False:
          for t in transcripts:
            axes[i].axvline(x = t[0], ymin = 0.5 if t[2] else 0, ymax = 1 if t[2] else 0.5, color = 'lawngreen')
            axes[i].axvline(x = t[1], ymin = 0.5 if t[2] else 0, ymax = 1 if t[2] else 0.5, color = 'lightcoral')


    #Plot whatever predictions are in the prediction list (designed to represent ML model predictions)
    if prediction:
      for i,p in enumerate(prediction):
        if p == 1:
          axes[i].axvline(x = i+1, color = 'green', alpha = 0.9, ymax = 0.3)
        elif p == 2:
          axes[i].axvline(x = i+1, color = 'red', alpha = 0.9, ymax = 0.3)
    
    #Plot whatever positions are in the peaks list
    if peaks:
      if type(peaks) != list: raise TypeError('Peaks must be a list of indices')
      for ind,j in enumerate(peaks):
        if j < 20: j = 25
        if peak_names:
          axes[i].axvline(x = j, color = peaks_color, alpha = 0.9, label = peak_names[ind])
        else:
          axes[i].axvline(x = j, color = peaks_color, alpha = peaks_alpha, ymax = 1)
    
    #Plot tis_term from the tis_term dataframe created by the proseq peak caller program
    if type(tis_term_df) == pd.DataFrame:
      if strand in [True, False]:
        tis_term_df = tis_term_df.loc[tis_term_df.Strand == ('Heavy' if strand else 'Light'), :]
      for _, row in tis_term_df.iterrows():
        checker = False
        loc = row.Position
        if loc < 20:
          loc += 20
          checker = True
        colour = 'green' if row.Type == 'TIS' else 'red'
        alpha = row.Confidence_score
        if not strand and transcript_overlap_mode:
          axes[i].axvline(x = loc, color = colour, alpha = alpha, ymax = 0.5)
        else:
          axes[i].axvline(x = loc, color = colour, alpha = alpha, ymax = 1)
        if peaks_annotate: # If peak annotation option is turned on, adds position values (x-axis) to the TIS/TERM.
          if plot_range[0] <= loc <= plot_range[1]:
            axes[i].text(x = loc + ((abs(plot_range[1] - plot_range[0])) * .005), y = y.max() * .7, s = loc - 20 if checker else loc, fontsize = 10, color = 'black', alpha = .5, fontstyle = 'oblique')
  
  plotsize = plot_range[1] - plot_range[0]
  if plotsize >= 1000:
    graphic_record = MitoTranslator().translate_record(record)
  else:
    graphic_record = MitoTranslatorT().translate_record(record)
  try:
    graphic_record = graphic_record.crop(plot_range)
    graphic_record.plot(ax = axes[0], with_ruler = False, annotate_inline = True, strand_in_label_threshold = 7)
  except AttributeError:pass
  if plotsize <= 120:
    graphic_record.plot_sequence(ax = axes[0], y_offset = 2)
  if log: axes[1].set_ylabel('log(Coverage)')
  plt.yscale(scale)
  if title_org:
    plt.suptitle(f'{org}{"_" if custom_title != "" else ""}{custom_title}', fontsize = 12)
  else:
    plt.suptitle(f'{custom_title}', fontsize = 12)
  axes[1].set_title(('Forward' if strand else 'Reverse') + ' strand')
  #axes[1].set_ylim(bottom = 0) #TODO Make sure this is a good idea
  if y_up_lim:
    axes[1].set_ylim(top = y_up_lim)
  plt.tight_layout()
  # Move legend outside of plot
  if legend_outside:
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
  else:
    axes[1].legend()
  # Move legend outside

  if return_ax:
    return fig, axes
  if savefig != '':
    fig.savefig(fname = savefig, dpi = 300, bbox_inches = 'tight')

def plot_gorder(org = '', ID = '', peaks = None, names = None, org_name = True, savefig = False):
  """
  Simple function to plot only the organism - without any data
  """

  if ID == '':
    check = record_check(org)
  else: 
    check = record_check(ID, mode = 'ID')
  if check == False:
    print(f'Organism {org} not found in the database')
    return
  record = SeqIO.read(f'{PATH}{slash}genbank_DB{slash}{org if org != "" else ID}.gbk', format = 'genbank')
  graphic_record = MitoTranslator().translate_record(record)
  try:
    title = record.annotations['references'][0].title + '\n' + record.annotations['references'][0].authors
  except(KeyError, IndexError):
    title = ''
  _, ax = plt.subplots(figsize = (15, 2.5))
  graphic_record.plot(ax = ax, with_ruler = True, annotate_inline = True, strand_in_label_threshold = 10)
  org = record.annotations['organism']
  if names:
    x_list = []
    for i, j in enumerate(peaks):
      # Add a dodge to the x-axis for each peak
      y = -3
      ax.axvline(x = j, color = 'black', alpha = 1, ymax = .2, linewidth = 1.9)
      ax.text(x = j, y = y, s = names[i], fontsize = 12, color = 'black', alpha = 1, fontstyle = 'normal', backgroundcolor = 'lightgrey')
      x_list.append(j)
  if org_name:
    plt.title(org, fontsize = 12, fontstyle = 'italic')
    if title != '':
      plt.suptitle(title, fontsize = 8)
  if savefig:
    plt.savefig(fname = os.path.join(PATH, 'figures', org + '.jpg'), dpi = 300)
  plt.tight_layout()

def plot_around_focal_point(ctrl, treat, tss_pos, left_region_size, right_region_size, pos_col = 'Position', cov_col = 'coverage', pos_prefix = 'pos_', neg_prefix = 'neg_', figsize = (14,4), treat_label = 'Hypoxia', ctrl_label = 'Normoxia', treat_color = 'red', ctrl_color = 'blue', plt_style = 'default', xlab = 'Position Relative to TSS', ylab = 'Coverage', savefig = False, fill = False, alpha = .7):
  """
  """
  plt.style.use(plt_style)

  pos_cov_col = pos_prefix + cov_col
  neg_cov_col = neg_prefix + cov_col
  mtDNA_len = int(ctrl.Position.max())
  left_pos = utils.mtdna_region(tss_pos, left_region_size, left = True, total = mtDNA_len, inclusive = False)
  right_pos = utils.mtdna_region(tss_pos, right_region_size, left = False, total = mtDNA_len, inclusive= False)

  # Filter df to only include the region of interest
  treat_reg_df = treat[treat.Position.isin(left_pos + right_pos)].copy()
  control_reg_df = ctrl[ctrl.Position.isin(left_pos + right_pos)].copy()

  # Set the TSS_pos as 0, the left_pos as negative values, and the right_pos as positive values
  treat_reg_df.loc[treat_reg_df.Position == tss_pos, 'Position'] = 0
  treat_reg_df.loc[treat_reg_df.Position.isin(left_pos), 'Position'] = range(-len(left_pos), 0)
  treat_reg_df.loc[treat_reg_df.Position.isin(right_pos), 'Position'] = range(1, len(right_pos) + 1)

  control_reg_df.loc[control_reg_df.Position == tss_pos, 'Position'] = 0
  control_reg_df.loc[control_reg_df.Position.isin(left_pos), 'Position'] = range(-len(left_pos), 0)
  control_reg_df.loc[control_reg_df.Position.isin(right_pos), 'Position'] = range(1, len(right_pos) + 1)
  
  # Convert negative strand coverage to negative values
  treat_reg_df = make_cols_neg(treat_reg_df)
  control_reg_df = make_cols_neg(control_reg_df)

  # Plot the data
  _, ax = plt.subplots(figsize = figsize)

  sns.lineplot(data = treat_reg_df, x = pos_col, y = pos_cov_col, ax = ax, label = treat_label, color = treat_color)
  sns.lineplot(data = treat_reg_df, x = pos_col, y = neg_cov_col, ax = ax, color = treat_color)
  sns.lineplot(data = control_reg_df, x = pos_col, y = pos_cov_col, ax = ax, label = ctrl_label, color = ctrl_color)
  sns.lineplot(data = control_reg_df, x = pos_col, y = neg_cov_col, ax = ax, color = ctrl_color)

  ax.set_xlabel(xlab)
  ax.set_ylabel(ylab)
  ax.legend()
  sns.despine()
  plt.tight_layout()

  if savefig:
    plt.savefig(fname = os.path.join(PATH, 'figures', f'{treat_label}_vs_{ctrl_label}.jpg'), dpi = 300)

