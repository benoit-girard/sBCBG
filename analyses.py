import pickle
import numpy as np

def main(path):
	NUCLEI = ['STN', 'GPe']
	file = open(path+'.txt', 'r')
	METRICS = pickle.load(file)
	print METRICS[0]['mat'].shape # (2, 6, 6, 15, 10, 1, 1)

	# Setting the positions and width for the bars
	pos = range(len(df['pre_score'])) 
	width = 0.25 
	    
	# Plotting the bars
	fig, ax = plt.subplots(figsize=(10,5))

	# Create a bar with pre_score data,
	# in position pos,
	plt.bar(pos, 
	        #using df['pre_score'] data,
	        df['pre_score'], 
	        # of width
	        width, 
	        # with alpha 0.5
	        alpha=0.5, 
	        # with color
	        color='#EE3224', 
	        # with label the first value in first_name
	        label=df['first_name'][0]) 

	# Create a bar with mid_score data,
	# in position pos + some width buffer,
	plt.bar([p + width for p in pos], 
	        #using df['mid_score'] data,
	        df['mid_score'],
	        # of width
	        width, 
	        # with alpha 0.5
	        alpha=0.5, 
	        # with color
	        color='#F78F1E', 
	        # with label the second value in first_name
	        label=df['first_name'][1]) 

	# Create a bar with post_score data,
	# in position pos + some width buffer,
	plt.bar([p + width*2 for p in pos], 
	        #using df['post_score'] data,
	        df['post_score'], 
	        # of width
	        width, 
	        # with alpha 0.5
	        alpha=0.5, 
	        # with color
	        color='#FFC222', 
	        # with label the third value in first_name
	        label=df['first_name'][2]) 

	# Set the y axis label
	ax.set_ylabel('Score')

	# Set the chart's title
	ax.set_title('Test Subject Scores')

	# Set the position of the x ticks
	ax.set_xticks([p + 1.5 * width for p in pos])

	# Set the labels for the x ticks
	ax.set_xticklabels(df['first_name'])

	# Setting the x-axis and y-axis limits
	plt.xlim(min(pos)-width, max(pos)+width*4)
	plt.ylim([0, max(df['pre_score'] + df['mid_score'] + df['post_score'])] )

	# Adding the legend and showing the plot
	plt.legend(['Pre Score', 'Mid Score', 'Post Score'], loc='upper left')
	plt.grid()
	plt.show()

		



def sortByDistance(matrices):
	reference = np.rot90([[1 if x+y > 5 else 0. for y in range(6)] for x in range(6)], k=3)
	reference = fil.gaussian_filter(reference, sigma=1.2) #apply gaussian filter
	reference = normalize(reference) #normalize
	[dic.update({'distance':np.sum(np.square(reference-dic['mat']))}) for dic in matrices]
	return sorted(matrices, key=lambda x: x['distance'])



def distFromRef(matrix):
	reference = np.rot90([[1 if x+y > 5 else 0. for y in range(6)] for x in range(6)], k=3)
	reference = fil.gaussian_filter(reference, sigma=1.2) #apply gaussian filter
	reference = normalize(reference) #normalize
	return np.sum(np.square(reference-matrix))



def normalize(matrix):
	return (matrix-np.min(matrix))/(np.max(matrix)-np.min(matrix))



if __name__ == '__main__':
	main('Sequence_5400')