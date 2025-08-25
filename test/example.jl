# Usage example:

resampler = BilinearBase(source_geo_def, target_geo_def, radius_of_influence)
get_bil_info!(resampler)

result = get_sample_from_bil_info(resampler, data)
