"""Recall useage: >python edit_filename.py"""

import os
for filename in os.listdir('.'):
	os.rename(filename, filename.replace('crazyname_', ''))
print('Done!')
