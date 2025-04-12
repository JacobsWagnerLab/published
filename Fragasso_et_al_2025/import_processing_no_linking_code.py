# -*- coding: utf-8 -*-
"""
Code for importing and analyse Omnipose masks and related phase and fluroescence images, after drift-correction by SuperSegger, with no linking.

Output: 
- Dictionary with cropped images for each cell, at each frame, for all channels (phase, fluor1, fluor2, ...) and corresponding masks
- Dataframe of cell features, e.g. area, time, mean fluorescences, correlations (SCF). Only 'static' features are estimated here based on the single frame.
    
@author: alessio fragasso
Tue Mar 25 18:02:23 2025
"""
import numpy as np
import pandas as pd
from skimage.measure import regionprops, find_contours
from scipy.stats import spearmanr
from scipy import ndimage
import matplotlib.pyplot as plt

"""
Extracts cell features from labelled masks for each frame.

Returns:
- DataFrame with cell_id, frame number, time_min, centroid, bounding box, area, major and minor axes.
"""
def extract_cell_data(mask, phase, fluor1, fluor2, cell_ch_dict, frame_number, frame_interval, cell_id_counter, 
                      pos, exp_name, rep, pad = 5, show=False, offs_x = 20, offs_y = 20):
    data = []
    props = regionprops(mask)
    for prop in props:
        centroid = prop.centroid
        minr, minc, maxr, maxc = prop.bbox
        x1, x2, y1, y2 = minc - pad, maxc + pad, minr - pad, maxr+ pad
        img_y, img_x = mask.shape
        if x1 < offs_x or y1 < offs_y or x2 > img_x-offs_x or y2 > img_y-offs_y:
            continue
        major_axis = prop.major_axis_length
        minor_axis = prop.minor_axis_length
        area = prop.area
        
        orientation_rad = props[0].orientation
        angle_deg = np.degrees(orientation_rad)

        angle_deg = np.abs(angle_deg)
        
        if area<500 or area > 2000:
            continue
        
        perimeter = prop.perimeter
        roundness = (4 * np.pi * area) / (perimeter**2)
        
        cell_mask = np.zeros(mask.shape)
        cell_mask[mask_temp==prop.label] = prop.label
        
        cropped_mask =cell_mask[y1:y2,x1:x2]
        cropped_phase = phase[y1:y2,x1:x2]
        cropped_fluor1 = fluor1[y1:y2,x1:x2]
        cropped_fluor2 = fluor2[y1:y2,x1:x2]

        if np.size(cropped_mask)!=np.size(cropped_fluor1):
            continue
        mean_fluor1 = np.mean(cropped_fluor1[cropped_mask>0])
        mean_phase = np.mean(cropped_phase[cropped_mask>0])
        sd_phase = np.std(cropped_phase[cropped_mask>0])
        mean_fluor2 = np.mean(cropped_fluor2[cropped_mask>0])

        mask_ero1 = ndimage.binary_erosion(cropped_mask, iterations=1)   #Erode the cell mask by px_erode pixels and take all pixels within that region
        mask_peri = np.logical_not(ndimage.binary_erosion(mask_ero1, iterations=2)) * mask_ero1
        phase_peri = cropped_phase[mask_peri>0].ravel()
        mean_phase_peri = np.mean(phase_peri)
        mean_phase_peri_high = np.mean(phase_peri[phase_peri>np.percentile(phase_peri,40)])
        
        if mean_phase_peri_high>3300:
            continue
        
        cv_phase_peri = np.std(cropped_phase[mask_peri>0].ravel())/mean_phase_peri 
        
        mask_eroded  = mask_ero1
        fluor1_roi = cropped_fluor1[mask_eroded > 0].ravel()
        fluor2_roi = cropped_fluor2[mask_eroded > 0].ravel()
        
        percentile_threshold = 50
        percentile_threshold_2 = 90
        
        # Calculate threshold values for each channel based on the chosen percentile
        threshold1 = np.percentile(fluor1_roi, percentile_threshold)
        threshold2 = np.percentile(fluor2_roi, percentile_threshold)
    
        fluor1_norm = fluor1_roi/np.max(fluor1_roi) 
        fluor2_norm = fluor2_roi/np.max(fluor2_roi)
        threshold3 = np.percentile(fluor1_norm+fluor2_norm, percentile_threshold_2)
        
        
        # Select pixels that are above the threshold in both channels
        valid_idx = (fluor1_roi > threshold1) & (fluor2_roi > threshold2)
        fluor1_valid = fluor1_roi[valid_idx]
        fluor2_valid = fluor2_roi[valid_idx]
        
        valid_idx_2 = (fluor1_norm+fluor2_norm)> threshold3
        fluor1_valid_2 = fluor1_roi[valid_idx_2]
        fluor2_valid_2 = fluor2_roi[valid_idx_2]
        
        
        
        # fluor1_roi = np.percentile(fluor1_roi,
        
        
        SCF_ori, _ = spearmanr(fluor1_roi, fluor2_roi)
        SCF, p_value = spearmanr(fluor1_valid, fluor2_valid)
        
        SCF_2, p_value = spearmanr(fluor1_valid_2, fluor2_valid_2)
     
        label_id = prop.label
        label_id = f"cell{label_id:07d}"
        cell_id_counter.append(f"cell{len(cell_id_counter) + 1:07d}")
        if label_id in cell_id_counter:
            label_id = cell_id_counter[-1]

        cell_id = label_id + '_'+pos + '_' +exp_name
        
        
        data.append({
            "cell_id": cell_id,
            'exp':exp_name,
            'pos':pos,
            'rep':rep,
            "frame": frame_number,
            "time_min": frame_number * frame_interval,
            "centroid_x": centroid[1],
            "centroid_y": centroid[0],
            'roundness': roundness,
            "bbox_x_min": x1,
            "bbox_y_min": y1,
            "bbox_x_max": x2,
            "bbox_y_max": y2,
            "r_off": np.array([x1+1,y1+1]),
            'rcm': np.array([centroid[1],centroid[0]]),
            'angle': angle_deg,
            "area": area,
            'area_um': area*(px_size**2),
            "major_axis_length": major_axis,
            "minor_axis_length": minor_axis,
            'mean_phase': mean_phase,
            'sd_phase': sd_phase,
            'mean_fluor1': mean_fluor1,
            'mean_fluor2': mean_fluor2,
            'SCF_ori': SCF_ori,
            'SCF': SCF
        })
        
        
        cell_ch_dict[cell_id] = cropped_mask, cropped_phase, cropped_fluor1, cropped_fluor2 
        
        
        if show:
            contours_ero = measure.find_contours(mask_eroded,level=0.5)
            contours = measure.find_contours(cropped_mask,level=0.5)
            
            fig = plt.figure(figsize=(25,10))
       
            fig.suptitle(cell_id + ', frame = ' + str(k))
            ax1 = fig.add_subplot(2,4,1)
            plt.title('mean phase = ' + str(round(mean_phase,1)) )
            plt.imshow(cropped_phase)
            for contour in contours:
                plt.plot(contour[:,1],contour[:,0], linewidth = 1, color = 'yellow',linestyle='--')
            
            ax1 = fig.add_subplot(2,4,2)
            plt.title('Area = ' + str(area) + ', roundness = ' + str(round(roundness,2)))
            plt.imshow(cropped_mask)
            
            ax1 = fig.add_subplot(2,4,3)
            plt.title('mean_phase_peri_high =' + str(round(mean_phase_peri_high ,2)))
            plt.imshow(mask_peri)
            
            ax1 = fig.add_subplot(2,4,4)
            plt.title('mean fluor1 = ' + str(round(mean_fluor1,1)))
            plt.imshow(cropped_fluor1)
            for contour in contours_ero:
                plt.plot(contour[:,1],contour[:,0], linewidth = 1, color = 'yellow',linestyle='--')
            
            ax1 = fig.add_subplot(2,4,5)
            plt.title('mean fluor2 = ' + str(round(mean_fluor2,1)))
            
            plt.imshow(cropped_fluor2*cropped_mask)
            for contour in contours_ero:
                plt.plot(contour[:,1],contour[:,0], linewidth = 1, color = 'yellow',linestyle='--')
            
            ax1 = fig.add_subplot(2,4,6)
            plt.title('SCF_ori = ' + str(round(SCF_ori,1)))
            plt.scatter(fluor1_roi, fluor2_roi, alpha = 0.5)
            plt.scatter(fluor1_valid, fluor2_valid, c='orange', alpha = 0.5)
            plt.scatter(fluor1_valid_2, fluor2_valid_2, c = 'green', alpha = 0.5)
            plt.xlabel('Fluor1')
            plt.ylabel('Fluor2')
            
            
            ax1 = fig.add_subplot(2,4,7)
            plt.title('SCF = ' + str(round(SCF,1)))
            plt.scatter(fluor1_valid, fluor2_valid)
            plt.xlabel('Fluor1')
            plt.ylabel('Fluor2')
            
            ax1 = fig.add_subplot(2,4,8)
            plt.title('SCF_2 = ' + str(round(SCF_2,1)))
            plt.scatter(fluor1_valid_2, fluor2_valid_2)
            plt.xlabel('Fluor1')
            plt.ylabel('Fluor2')
            plt.tight_layout()
            plt.show()

    df = pd.DataFrame(data)
    return df, cell_ch_dict, cell_id_counter    



'''Start of the code that runs through the folders'''

'''
Input global variables
'''
fr = 1    # [min] time between frame
px_size = 0.065841  # um/px   (pixel size in alpha scope, 100x)
px_erode = 4                # number of pixels to erode in the cell mask to identify ROI for correlation
offs_x = 10                # distance from image edges to exclude from analysis
offs_y = 10
save = True
bkg_sub = True
show = False

pos_list = create_pos_list(1, 20, double_digit = True) 
ch_list = ['phase', 'fluor1', 'fluor2','fluor3']
ch_fluor_list = [s for s in ch_list if 'fluor' in s]

'''Specify paths'''
master_folder_path = #r"path_to_strain_master_folder"   # from mothership2

'''Insert list of AMP (conditions) to process'''
AMP_list = [ 
    ## AMP_1,
    ## AMP_2,
    ## AMP_3,
    ## ... 
    ]

'''Insert list of lists of frames at which AMP was injected. One list for each AMP condition that contains fr_inj for each replicate'''

fr_inj_master_list = [
   ##[fr_inj_rep_1, fr_inj_rep_2, fr_inj_rep_3],   # for AMP_1  
   ##[fr_inj_rep_1, fr_inj_rep_2, fr_inj_rep_3],   # for AMP_2  
   ##[fr_inj_rep_1, fr_inj_rep_2, fr_inj_rep_3],   # for AMP_3
    ]


for l in range(len(AMP_list)):
    AMP_name = AMP_list[l]
    fr_inj_list = fr_inj_master_list[l]
    rep_list = [s for s in os.listdir(master_folder_path+'/'+AMP_name) if '.DS_Store' not in s]
    # rep_list = ['250201_4_5per']
    
    for s in range(len(rep_list)):
        # s=2
        rep_name = rep_list[s]
        exp_name = AMP_name + '_' + rep_name + '_' + os.path.basename(master_folder_path)
        experiment_path = master_folder_path+'/'+AMP_name+'/'+rep_name
        print(exp_name)
        if len(fr_inj_list)==0:
            fr_inj = 0
        else:
            fr_inj = fr_inj_list[s]
        for pos in pos_list:
            '''
            Get paths to all folders
            '''
            
            folder_path = experiment_path + '/' +\
                pos  # path to position
            cells_path = folder_path+"/cell"
            mask_path = folder_path + '/masks'
            
            pos_num = re.findall(r'\d+', pos)[0]
         
            if os.path.exists(mask_path):
                masks_list = os.listdir(mask_path)
            else: continue
            print(pos)
            paths_dict = {}
            file_list_dict = {}
            
            for ch in ch_list:   
                paths_dict[ch] = folder_path + '/'+ ch
                file_list_dict[ch] = os.listdir(paths_dict[ch])
                    
            phase_list =file_list_dict['phase']
            phase_path = paths_dict['phase']
            '''
            Initialize dictionaries and lists
            '''
            cell_stack_dict = {}  # dictionary of cells parameters over time
            cell_ch_dict = {}  # dictionary of cells phase and mask over time
            cell_features_all_df = pd.DataFrame()   
            cell_id_counter = []  # Store unique cell IDs across frames
            '''Subtract background from fluorescence images and store in temporary dictionaries'''
            fluor_bkg_sub_dict = {}
            fluor_mask_dict = {}
            mask_dict = {}
            phase_dict = {}
            img_size = {}

            print('..Subtract background and extract cell features at pos ' + pos + ' for experiment '+ exp_name)
            fluor1_list = file_list_dict['fluor1']
            fluor1_path = paths_dict['fluor1']
            fluor2_list = file_list_dict['fluor2']
            fluor2_path = paths_dict['fluor2']
            for k in range(len(masks_list)):
                if AMP_name == '1_5T-CecA_1uM':
                    if k>14:
                        continue
                if k%20==0: 
                    print(k)
                mask_temp = iio.imread(mask_path+'/'+masks_list[k])
                phase_temp = iio.imread(phase_path+'/'+phase_list[k])
                fluor1_temp = iio.imread(fluor1_path + '/'+fluor1_list[k])
                fluor2_temp = iio.imread(fluor2_path + '/'+fluor2_list[k])
            
                bs = int(fluor1_temp.shape[0]/16)
                fluor1_bkg_sub, _, _ = local_bkg_sub_cp(
                    fluor1_temp, mask_temp, pos, exp_name, frame=str(k), box_size=bs, dilations=15, sigma_=60, show=False)
                fluor2_bkg_sub, _, _ = local_bkg_sub_cp(
                    fluor2_temp, mask_temp, pos, exp_name, frame=str(k), box_size=bs, dilations=15, sigma_=60, show=False)

                df_temp, cell_ch_dict, cell_id_counter = extract_cell_data(mask_temp, phase_temp, fluor1_bkg_sub, fluor2_bkg_sub, cell_ch_dict, k, fr,
                                                                           cell_id_counter, pos, exp_name,rep_name, pad = 5, show=show)
                
                df_temp['fr_inj'] = fr_inj
                cell_features_all_df = pd.concat([cell_features_all_df,df_temp])
                      
            if save and len(cell_features_all_df)>0:
                cell_ch_dict_name = pos + '_' + exp_name + '_cell_ch_bkg_sub'
                df_basename = pos + '_' + exp_name +'_cell_features_df'
                save_dict_path = experiment_path + '/output'
                if not os.path.exists(save_dict_path):
                    os.makedirs(save_dict_path)
                dict_cell_path = os.path.join(save_dict_path, cell_ch_dict_name+'.pkl')
                with open(dict_cell_path, 'wb') as pickle_file:
                    pickle.dump(cell_ch_dict, pickle_file)
                cell_features_all_df.to_pickle(save_dict_path+'/'+df_basename+'.pkl')
                count += len(cell_features_all_df['cell_id'].unique())
                print('At position '+pos+' there are ' +
                      str(len(cell_features_all_df['cell_id'].unique()))+' cells')
                del cell_ch_dict
                gc.collect()
    print('In total there are '+
          str(count)+' cells')  
            
    

