import numpy as np

def crop_array(samples, length):
    """
    Array will have a exact multiple length of length

    Parameters:
    - samples: array of the time samples
    - length: desired block length which will be uses as a reference for
    cropping

    Return: 
    - cropped samples
    - number of chunks
    """
    size = np.prod(samples.shape)
    # make 1-d array
    samples = samples.reshape([size, 1])
    # determine the multiple of length
    num_chunks = int(size / length)
    # return copped version and the number of chunks
    return samples[0:num_chunks*length], num_chunks

def block_averaging(samples, block_length):
    """
    Return a list of the averages of block_length made from samples

    Parameters:
    - samples: array of the time samples
    - block_length: length of the blocks used for averaging

    Return:
    - averages of the blocks
    - number of blocks
    """
    # crop to fit the desired length for making blocks
    samples, num_blocks = crop_array(samples, block_length)
    # each block is one row
    samples = samples.reshape([num_blocks, block_length])
    # compute the sums over the rows
    block_averages = samples.mean(axis=1)
    return block_averages, num_blocks

def block_variance(samples, block_length, time_average, stdout=True):
    """
    Return the variance of the block averages

    Parameters:
    - samples: array of the time samples
    - block_length: length of the blocks used for averaging
    - time_average: reference mean value of the samples
    - stdout: flag to disable the printing of debug information

    Return:
    - variance of the block averages
    - number of blocks
    """
    # compute the block averages
    block_averages, num_blocks = block_averaging(samples, block_length)
    # print info if wanted
    if stdout:
        print(f"block length: {block_length} \t, num blocks: {num_blocks} \t, used samples: {num_blocks * block_length}")
    # compute the standard error of the block averages
    block_var = np.sum((block_averages - time_average)**2) / num_blocks

    return block_var, num_blocks


def time_estimate(samples, block_length):
    """
    Give the time average and the standard error

    Parameters:
    - samples: array of the time samples
    - block_length: length of the blocks to be used for computing the standard
    error
    """
    time_average = np.mean(samples)
    block_var = block_variance(samples, block_length, time_average, stdout=False)
    std_err = np.sqrt(block_var[0] / (block_var[1] - 1))
    print(f"{time_average:.3f} +/- {std_err:.3f}")
    return time_average, std_err