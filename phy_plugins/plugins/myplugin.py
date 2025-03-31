from phy import IPlugin
from phy.cluster.views import ManualClusteringView  # Base class for phy views
from phy.plot.plot import PlotCanvasMpl  # matplotlib canvas
from numpy import genfromtxt
import matplotlib
import matplotlib.pyplot as plt
import colorcet as cc
import numpy as np
import warnings
import time
import sys
import os
from phy.utils.color import selected_cluster_color

import scipy.io

warnings.filterwarnings("ignore")

axis_list = []
plot_handles = []
nodata=False
growth = 1
decay = 20
A = np.empty((0, 2))

for i in range(101):  # Equivalent to 0:1:100 in MATLAB
    A = np.vstack((A, [i, (1 - np.exp(i * -1 / growth)) * (np.exp(i * -1 / decay))]))
    B = np.cumsum(A, axis=0)
    C = B[-1, -1]
    KERNEL = np.column_stack((A[:, 0], A[:, 1] / C))
    
def _make_default_colormap():
    """Return the default colormap, with custom first colors."""
    colormap = np.array(cc.glasbey_bw_minc_20_minl_30)
    # Reorder first colors.
    colormap[[0, 1, 2, 3, 4, 5]] = colormap[[3, 0, 4, 5, 2, 1]]
    # Replace first two colors.
    colormap[0] = [0.03137, 0.5725, 0.9882]
    colormap[1] = [1.0000, 0.0078, 0.0078]
    return colormap


class RasterPSTH_class(ManualClusteringView):
    plot_canvas_class = PlotCanvasMpl  # use matplotlib instead of OpenGL (the default)
    def __init__(self, c=None):
        """features is a function (cluster_id => Bunch(data, ...)) where data is a 3D array."""
        super(RasterPSTH_class, self).__init__()
        self.controller = c
        self.model = c.model
        self.supervisor = c.supervisor
        self.cmap=_make_default_colormap()

        self.event_idx = 0 # Determines the alignment, target onset, saccade onset, arm-movement onset
        
        # Initialize an array for numbers and a dictionary for descriptions
        conds = []
        conds_des = {}
        
        print(os.path.dirname(c.model.dir_path))
        print(c.model.dir_path)
        print(os.path.join(os.path.dirname(c.model.dir_path), 'condition.txt'))
              
        
        mat = scipy.io.loadmat(os.path.dirname(c.model.dir_path) + '\\conditions.mat')
        self.conds = mat['conditions'][0]  # This will give you the data as a numpy array
        # # Open and read the text file
        # with open(os.path.join(os.path.dirname(c.model.dir_path), 'condition.txt'), 'r') as file:
        #     print('hi')
        #     for line in file:
        #         # Split the line by tab character
        #         parts = line.strip().split('\t')
                
        #         # Store the number and description
        #         number = int(parts[0])
        #         description = parts[1]
                
        #         # Append the number to the list and add to the dictionary
        #         conds.append(number)
        #         conds_des[number] = description
        
        mat = scipy.io.loadmat(os.path.dirname(c.model.dir_path) + '\\event_times_r.mat')
        self.events_r = mat['event_times_r']  # This will give you the data as a numpy array

        mat = scipy.io.loadmat(os.path.dirname(c.model.dir_path) + '\\event_times_l.mat')
        self.events_l = mat['event_times_l']  # This will give you the data as a numpy array
            
        self.conds_des = conds_des
        # self.conds = conds
        self.cond_idx = 0
        
        # self.events_rr = events_r
        # self.events_ll = events_l
        
    def on_request_similar_clusters(self,cid=None):
        self.on_select()

    def on_select(self, cluster_ids=(), **kwargs):
        if nodata:
            return

        global axis_list,plot_handles
        cluster_ids=self.supervisor.selected
        self.cluster_ids=cluster_ids
        self.nclusts = len(cluster_ids)

        if axis_list:
            axis_diff = (self.nclusts-len(axis_list)//2)*2
            if axis_diff<0:
                axis_list = axis_list[0:(len(axis_list))+axis_diff]
                plot_handles = plot_handles[0:len(plot_handles)+axis_diff//2]

        # We don't display anything if no clusters are selected.
        if not cluster_ids:
            return

        for i,d in enumerate(np.arange(start=1,stop=self.nclusts*2+1)):
            if d%2 == 0:
                setattr(self,'canvas.ax'+str(d), plt.subplot(2*self.nclusts,1,d,sharex=axis_list[i-1]))
            else:
                setattr(self,'canvas.ax'+str(d), plt.subplot(2*self.nclusts,1,d))

            if (len(axis_list)-1)<i:
                axis_list.append(getattr(self,'canvas.ax'+str(d)))
            else:
                axis_list[i]=(getattr(self,'canvas.ax'+str(d)))
            axis_list[i].cla()

        t1=time.time()
        ttime=0
        
        print(self.event_idx)
        for i,d in enumerate(cluster_ids):
            flag = 0
            # _________________________________ right _______________________________________________________
            if(self.events_r.size > 0):
                events_r = self.events_r[0, self.conds[self.cond_idx]-1][self.event_idx]
                events_r = events_r.reshape(-1)
                spdn, rasters,activity,yrast,ntrials,nevents=self.get_spikes(d, events_r[0:-1])
                
                t=(time.time())
                axis_list[i*2+1].scatter(rasters,yrast,c=selected_cluster_color(i),
                    marker='|',s=np.ones(len(rasters))*1,alpha=0.8)
                axis_list[i*2+1].axvline(x=0,color='white',alpha=.5)
                #hist, bins = np.histogram(activity,weights=np.ones(nevents)*(50/ntrials),range=(-1,1),bins=100)
                # axis_list[i*2].plot(bins[:-1],hist,color=self.cmap[i])
                
                axis_list[i*2].plot(np.linspace(-0.5, 0.5, 1001), spdn,c=selected_cluster_color(i), label = 'Right')
                axis_list[i*2].legend(loc='upper left', bbox_to_anchor=(1, 1))
                axis_list[i*2].axvline(x=0,color='white',alpha=.5)
                axis_list[i*2].set_xticks(np.arange(-0.2, 0.5 + 0.05, 0.05))
                axis_list[i*2].set_xlim(left=-0.2,right=0.5)
                axis_list[i*2+1].set_ylim(bottom=0,top=ntrials)
                #axis_list[i*2].set_ylim(bottom=0,top=spM+spM/10)
                self.fix_axis(axis_list[i*2],10)
                self.fix_axis(axis_list[i*2+1],10)
                flag = 1
            
            # _________________________________ left _______________________________________________________
            if(self.events_l.size > 0):
            
                events_l = self.events_l[0, self.conds[self.cond_idx]-1][self.event_idx]
                events_l = events_l.reshape(-1)
                spdn2, rasters2,activity2,yrast2,ntrials2,nevents=self.get_spikes(d, events_l[0:-1])

                t=(time.time())
                #print(yrast2+np.max(yrast))
                if flag:
                    axis_list[i*2+1].scatter(rasters2,yrast2+np.max(yrast),c=selected_cluster_color(i+1),
                        marker='|',s=np.ones(len(rasters2))*1,alpha=0.8)
                else:
                    axis_list[i*2+1].scatter(rasters2,yrast2,c=selected_cluster_color(i+1),
                        marker='|',s=np.ones(len(rasters2))*1,alpha=0.8)

                axis_list[i*2+1].axvline(x=0,color='white',alpha=.5)
                #hist, bins = np.histogram(activity,weights=np.ones(nevents)*(50/ntrials),range=(-1,1),bins=100)
                # axis_list[i*2].plot(bins[:-1],hist,color=self.cmap[i])
                #axis_list[i*2].plot([0.25, 0.5, 0.75], [2, 3, 4])
                
                axis_list[i*2].plot(np.linspace(-0.5, 0.5, 1001), spdn2,c=selected_cluster_color(i+1), label = 'Left')
                axis_list[i*2].legend(loc='upper left', bbox_to_anchor=(1, 1))
                axis_list[i*2].axvline(x=0,color='white',alpha=.5)
                axis_list[i*2].set_xticks( np.arange(-0.2, 0.5 + 0.05, 0.05))#  np.arange(x, y + step, step)
                axis_list[i*2].set_xlim(left=-0.2,right=0.5)
                if flag:
                    axis_list[i*2+1].set_ylim(bottom=0,top=ntrials2+np.max(yrast))
                else:
                    axis_list[i*2+1].set_ylim(bottom=0,top=ntrials2)

                #axis_list[i*2].set_ylim(bottom=0,top=spM+spM/10)
                self.fix_axis(axis_list[i*2],10)
                self.fix_axis(axis_list[i*2+1],10)
        
        cond = self.conds[self.cond_idx]

        if cond == 1:
            txt = 'EYE_EyeOnly'
        elif cond == 2:
            txt = 'HAND_Hand'
        elif cond == 3:
            txt = 'EYEHAND_EyeOnly'
        elif cond == 4:
            txt = 'EYEHAND_EyeHand'
        elif cond == 5:
            txt = 'ETP: Single Target'
        elif cond == 6:
            txt = 'ETP: Eye Only'
        elif cond == 7:
            txt = 'Mapping'
        elif cond == 8:
            txt = 'Delayed Task'
        elif cond == 9:
            txt = 'ETP: Single Target, Static'
        else:
            txt = ''
                
            
        if(self.event_idx == 0):
            axis_list[0].set_title(txt + " - Aligned to Target Onset")
        elif(self.event_idx == 1):
            axis_list[0].set_title(txt + " - Aligned to Saccade Onset")
        elif(self.event_idx == 2):
            axis_list[0].set_title(txt + " - Aligned to Arm Onset")

        # print(ttime,'plotting')
        t=time.time()
        # Use this to update the matplotlib figure.
        self.canvas.show()
        
        # print(time.time()-t,'update')
        # print(time.time()-t1,'total')
        return


    def fix_axis(self,ax,textsize):
        ax.tick_params(axis='x', labelsize=textsize)
        ax.tick_params(axis='y', labelsize=textsize)
        ax.xaxis.label.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.grid(False)
    
    def get_spikes(self,clust, events):
        spikes = self.model.get_cluster_spikes(clust)
        spike_times =np.array(self.model.spike_times[spikes])
        rasters = np.array([])
        yrast = np.array([])
        activity = []
        last_ind = 0
        max_spike_times = np.amax(spike_times)
        min_spike_times = np.amin(spike_times)
        spden_aligned_mat = []
        i2 = 0
        for i,d in enumerate(events):
            if d<max_spike_times and d>min_spike_times:
                st = spike_times[last_ind:]
                temp1 = st-(d+0.5)
                ind1 = np.abs(temp1).argmin()
                if temp1[ind1]>0:
                    ind1-=1
                temp2 = st-(d-0.5)
                ind2 = np.abs(temp2).argmin()
                if temp2[ind2]<0:
                    ind2+=1
                temp=st[ind2:ind1]-d
                i2 = i2 + 1
                #print(i2, len(temp))
                #if(len(temp) > 0):
                
                last_ind=ind1
                rasters=np.append(rasters,temp)
                #print(len(temp))
                yrast=np.append(yrast,np.ones(len(temp))*i2)
                activity.extend(temp)
                #print(temp)
                # Assuming spike_times_aligned is a list or array of spike times
                CC = np.zeros(1001)  # Row of zeros for timespan of trial
                
                #print(np.round(1000*temp))
                #print(np.round(1000*temp).astype(int))
                CC[np.round(1000*(temp)).astype(int) + 500] = 1  # Replaces 0 with 1 at spike_times
                #print(np.sum(CC))   
                #print(CC.shape)
                #print(KERNEL.shape)
                spden_aligned = np.convolve(CC, KERNEL[:, 1], mode='full') * 1000
                
                #print(spden_aligned.shape)
                spden_aligned = spden_aligned[:1001]
                spden_aligned_mat.append(spden_aligned)
                
                #print(len(spden_aligned_mat))
                #print(spden_aligned_mat)
                #spden_aligned_mat[i, :] = spden_aligned
                
        spdn = np.sum(spden_aligned_mat, axis=0) / len(events)
        #spdn = np.mean(spden_aligned_mat, axis=0)      
        #print(spdn.shape)
        #print(spdn)
        #print(len(spdn))
        return spdn, rasters,np.array(activity),yrast,np.amax(yrast),len(activity)

from phy import IPlugin, connect
from phy import emit
events_r = 0
events_l = 0
class RasterPSTH(IPlugin):
    def attach_to_controller(self, controller):
        def create_PSTH_view():
                """A function that creates and returns a view."""
                view = RasterPSTH_class(controller)

                @connect(sender=view)
                def on_view_attached(view_, gui):
                    # NOTE: this callback function is called in PSTHView.attach().
                    # @view.actions.add(prompt=True, prompt_default=lambda: str(view.t_pre))
                    # def change_t_pre(t_pre):
                    #     """Change the t_pre in millisecond displayed in the PSTHView (enter positive number)."""
                    #     view.t_pre = -t_pre
                    #     view.on_select(view.cluster_ids) 

                    # @view.actions.add(prompt=True, prompt_default=lambda: str(view.t_post))
                    # def change_t_post(t_post):
                    #     """Change the t_post in millisecond displayed in the PSTHView."""
                    #     view.t_post = t_post
                    #     view.on_select(view.cluster_ids) 

                    # @view.actions.add(prompt=True, prompt_default=lambda: str(view.binwidth))
                    # def change_binwidth(binwidth):
                    #     """Change the binwidth in millisecond displayed in the PSTHView."""
                    #     view.binwidth = binwidth
                    #     view.on_select(view.cluster_ids) 

                    @view.actions.add(shortcut='e')
                    def next_event():
                        """Change to the next event displayed in the PSTHView."""
                        view.event_idx = np.mod(view.event_idx+1,view.events_r[0, view.conds[view.cond_idx]-1].shape[0])  
                        print(view.event_idx)
                        view.on_select() 
                        
                    @view.actions.add(shortcut='q')
                    def next_condition():
                        """Change to the next condition displayed in the PSTHView."""
                        view.cond_idx = np.mod(view.cond_idx+1,view.conds.shape[0])  
                        print(view.event_idx)
                        view.on_select() 


                return view

        controller.view_creator['RasterPSTH'] = create_PSTH_view 
                
                

   
# import from plugins/cluster_metadata.py
"""Show how to save the best channel of every cluster in a cluster_channel.tsv file when saving.

Note: this information is also automatically stored in `cluster_info.tsv` natively in phy,
along with all values found in the GUI cluster view.

"""

import logging

from phy import IPlugin, connect
from phylib.io.model import save_metadata

logger = logging.getLogger('phy')

class SaveMetadataPlugin(IPlugin):
    def attach_to_controller(self, controller):
        @connect
        def on_gui_ready(sender, gui):

            @connect(sender=gui)
            def on_request_save(sender):
                """This function is called whenever the Save action is triggered."""

                # We get the filename.
                filename = controller.model.dir_path / 'cluster_channel.tsv'

                # We get the list of all clusters.
                cluster_ids = controller.supervisor.clustering.cluster_ids

                # Field name used in the header of the TSV file.
                field_name = 'channels'

                # NOTE: cluster_XXX.tsv files are automatically loaded in phy, displayed
                # in the cluster view, and interpreted as cluster labels, *except* if their
                # name conflicts with an existing built-in column in the cluster view.
                # This is the case here, because there is a default channel column in phy.
                # Therefore, the TSV file is properly saved, but it is not displayed in the
                # cluster view as the information is already shown in the built-in channel column.
                # If you want this file to be loaded in the cluster view, just use another
                # name that is not already used, like 'best_channel'.

                # Dictionary mapping cluster_ids to the best channel id.
                metadata = {
                    cluster_id: controller.get_best_channel(cluster_id)
                    for cluster_id in cluster_ids}

                # Save the metadata file.
                save_metadata(filename, field_name, metadata)
                logger.info("Saved %s.", filename)
                