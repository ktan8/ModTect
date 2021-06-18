def calculate_median(values):
	values.sort()
	num_of_values = len(values)
	
	# odd case
	if num_of_values % 2 == 1:
		reqr_index = (num_of_values - 1) / 2
		reqr_value = values[reqr_index]

		median = reqr_value
		return median

	# Even case
	else:
		reqr_index_small = num_of_values / 2 - 1
		reqr_index_big = num_of_values / 2
		reqr_value_small = values[reqr_index_small]
		reqr_value_big = values[reqr_index_big]

		median = (float(reqr_value_small) + float(reqr_value_big)) / 2

		return median



def calculate_median_absolute_deviation(median_value, values):
	deviations = []
	for val in values:
		deviation_val = abs(float(val) - float(median_value))
		deviations.append(deviation_val)

	median_absolute_deviation = calculate_median(deviations)

	return median_absolute_deviation



# data = [20,30,31,23,24]
# data = [,30,31,23,24]
# medianVal = calculate_median(data)
# print calculate_median_absolute_deviation(medianVal, data)
