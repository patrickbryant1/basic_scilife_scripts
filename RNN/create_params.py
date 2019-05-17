#! /usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Prints files with different combinations of parameters to be used to define different models.
'''

num_res_blocks=[1,2,3,4,5]
base_epochs=[10,20,30]
finish_epochs=[3]
num_nodes=[128,192,256]
embedding_size=[5,10]
drop_rate=[0.5,0.7]
find_lr=[0]

for block in num_res_blocks:
	for base in base_epochs:
		for finish in finish_epochs:
			for nodes in num_nodes:
				for embedding in embedding_size:
					for drop in drop_rate:
						for lr in find_lr:
							name = str(block)+'_'+str(base)+'_'+str(finish)+'_'+str(nodes)+'_'+str(embedding)+'_'+str(drop)+'_'+str(lr)+'.params'
							with open(name, "w") as file: 
															
								file.write('num_res_blocks='+str(block)+'\n')
								file.write('base_epochs='+str(base)+'\n')
								file.write('finish_epochs='+str(finish)+'\n')
								file.write('num_nodes='+str(nodes)+'\n')
								file.write('embedding_size='+str(embedding)+'\n')
								file.write('drop_rate='+str(drop)+'\n')
								file.write('find_lr='+str(lr)+'\n')

