#! /usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Prints files with different combinations of parameters to be used to define different models.
'''

num_res_blocks=[1,2,3,4,5]
base_epochs=[10,20,40]
finish_epochs=[3]
filters = [100,300,600,900]
dilation_rate = [3,6]


for block in num_res_blocks:
	for base in base_epochs:
		for finish in finish_epochs:
			for filter in filters:
				for dr in dilation_rate:
					name = str(block)+'_'+str(base)+'_'+str(finish)+'_'+str(filter)+'_'+str(dr)+'.params'
					with open(name, "w") as file:
						file.write('num_res_blocks='+str(block)+'\n')
						file.write('base_epochs='+str(base)+'\n')
						file.write('finish_epochs='+str(finish)+'\n')
						file.write('filters='+str(filter)+'\n')
						file.write('dilation_rate='+str(dr)+'\n')
