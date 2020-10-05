import numpy as np
from copy import copy as copy


def get_pixel_shells(max_rad, min_rad=0):
    
    arr = np.arange(2*max_rad+1) - max_rad
    X,Y = np.meshgrid(arr, arr)
    R = np.round(np.sqrt(X**2 + Y**2)).astype(int)
    
    indices = {}
    for i in range(min_rad, max_rad+1):
        indices[i] = np.where(R == i)
        indices[i] = (indices[i][0]-max_rad, indices[i][1]-max_rad)
    
    return indices


def merge_data(images, nanMaps, centers, std_cut, maxPixel,
    rad_inds, ref_rad_pixs=None, ref_rad_nans=None, loop_std_cut=1):
    """
    Must have dimension order [NtimeSteps, Nscans, Nrows, Ncols]
    By specifying ref_rad_pixels it will subtract the reference
    for rad in rads | Could make a radial loop if needed in the future
    """
    
    merged_centerR = maxPixel - 1
    merged_centerC = maxPixel - 1
    merged_image = np.ones((len(images), 2*maxPixel-1, 2*maxPixel-1))*np.nan
    merged_image_std = np.ones_like(merged_image)*np.nan
    merged_image_counts = np.ones_like(merged_image)*np.nan
    for rad in range(maxPixel):
        pixels, nans, angles = [], [], []
        angles = np.arctan2(rad_inds[rad][1], rad_inds[rad][0])
        
        for itm in range(len(images)):
            scn_inds, row_inds, col_inds = [], [], []
            for isc in range(images[itm].shape[0]):
                scn_inds.append(np.ones(rad_inds[rad][0].shape[0], dtype=int)*isc)
                row_inds.append(rad_inds[rad][0]+centers[itm][isc,0])
                col_inds.append(rad_inds[rad][1]+centers[itm][isc,1])
            scn_inds = np.concatenate(scn_inds)
            row_inds = np.concatenate(row_inds)
            col_inds = np.concatenate(col_inds)
            img_pixels = np.reshape(copy(images[itm][scn_inds,row_inds,col_inds]), (images[itm].shape[0], -1))
            img_nans   = np.reshape(copy(nanMaps[itm][scn_inds,row_inds,col_inds]), (images[itm].shape[0], -1))
            
            if ref_rad_pixs is not None:
                img_pixels -= ref_rad_pixs[rad]
                img_nans   = (img_nans + np.expand_dims(ref_rad_nans[rad], 0)).astype(bool)
                images[itm][scn_inds,row_inds,col_inds] -= np.tile(ref_rad_pixs[rad], images[itm].shape[0])
            
            
            mask = np.ones_like(img_pixels)
            mask[img_nans.astype(bool)] = np.nan
            for istd in range(loop_std_cut):
              mean = np.nanmean(img_pixels*mask, axis=0)
              std  = np.nanstd(img_pixels*mask, axis=0)
              std[std==0.0] = np.sqrt(np.abs(mean[std==0.0]))
          
              # Cut pixels 
              mask[np.abs(mean-img_pixels) > std_cut*std] = np.nan
              img_nans[np.abs(mean-img_pixels) > std_cut*std] = 1
                
            
            # Merge
            merged_img_pixels = np.nanmean(img_pixels*mask, axis=0)
            merged_img_std    = np.nanstd(img_pixels*mask, axis=0)
            merged_img_counts = np.nansum(mask, axis=0)
            merged_img_std    /= np.sqrt(merged_img_counts)
            merged_img_std[np.isnan(merged_img_std)] = 0
            
            nanMaps[itm][scn_inds,row_inds,col_inds] = img_nans.flatten()
            merged_image[itm,merged_centerR+rad_inds[rad][0],merged_centerC+rad_inds[rad][1]] = copy(merged_img_pixels[:])
            merged_image_std[itm,merged_centerR+rad_inds[rad][0],merged_centerC+rad_inds[rad][1]] = copy(merged_img_std[:])
            merged_image_counts[itm,merged_centerR+rad_inds[rad][0],merged_centerC+rad_inds[rad][1]] = copy(merged_img_counts[:])
        
    
    return merged_image, merged_image_std, merged_image_counts, nanMaps, images
