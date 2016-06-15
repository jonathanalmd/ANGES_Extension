# ANGES 1.0, reconstruction ANcestral GEnomeS maps
# May 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import copy

#######################################################
#   pqtree.py
#
#    Classes to encode pq-trees.
#
#######################################################


# a PQCR tree
class PQTree:
	# creates a new PQTree
	# name - str: name of ancestor species the PQCR tree represents
	def __init__(self, name):
		self.name = name
		self.head = TreeNode(None)		# the head of the PQTree
	#enddef

	# returns a string of the PQTree in its linear representation
	# return - str
	def __str__(self):
		s = '>' + self.name + '\n'

		if self.head.value != 'P':
			s += '#CAR0\n'
			s += self.head.str_self() + '\n'
		else:
			car = self.head.child
			i = 1

			while car != None:
				s += '#CAR' + str(i) + '\n'
				s += car.str_self() + '\n'

				car = car.brother
				i += 1
			#endfor
		#endif

		return s
	#enddef

	def to_file_name(self, file_name):
		file_stream = open(file_name, 'w')

		file_stream.write('>' + self.name + '\n')

		if self.head.value != 'P':
			file_stream.write('#CAR0\n')
			self.head.write_self(file_stream)
			file_stream.write('\n')
		else:
			car = self.head.child
			i = 1

			while car != None:
				file_stream.write('#CAR' + str(i) + '\n')
				car.write_self(file_stream)
				file_stream.write('\n')

				car = car.brother
				i += 1
			#endfor
		#endif


		file_stream.close()
	#enddef

	def __copy__(self):
		tre = PQRTree()

		tre.head = copy.copy(self.head)

		return tre
	#enddef

	@staticmethod
	def from_file_name(file_name):
		f = open(file_name, 'r')

		tree = PQTree('')
		tree.head = TreeNode('P')
		node = tree.head

		for line in f:
			stripped = line.strip()

			if stripped.startswith('>'):
				tree.name = stripped[1:]
				continue

			if stripped == '' or stripped.startswith('#'):
				continue

			if node == tree.head:
				tree.head.child, string = TreeNode.from_string(line)
				node = tree.head.child
			else:
				node.brother, string = TreeNode.from_string(line)
				node.brother.last = node
				node = node.brother

		return tree


# node in a PQCR tree
class TreeNode:
	# creates a new TreeNode with value, val
	# val - anything
	def __init__(self, val):
		self.value = val		# the node's value
		self.child = None		# the node's child
		self.brother = None		# the node's brother
		self.last = None		# the previous brother (only needed
								# for joining partitions)
	#enddef

	# returns a string of the node and all its brothers and children
	# return - str
	def __str__(self):
		if self.child != None:
			# the string representation of the TreeNode
			s = '_' + str(self.value) + ' ' + self.child.__str__() + str(self.value) + '_ '
		else:
			s = str(self.value) + ' '
		#endif

		if self.brother != None:
			s += str(self.brother)
		#endif

		return s
	#enddef

	# returns a string of the node and all its children
	# return - str
	def str_self(self):
		if self.child != None:
			# the string representation of the TreeNode
			s = '_' + str(self.value) + ' ' + self.child.__str__() + str(self.value) + '_ '
		else:
			s = str(self.value)
		#endif

		return s
	#enddef

	def write(self, file_stream):
		if self.child != None:
			file_stream.write('_' + str(self.value) + ' ')
			self.child.write(file_stream)
			file_stream.write(str(self.value) + '_ ')
		else:
			file_stream.write(str(self.value) + ' ')
		#endif

		brother = self.brother

		while brother != None:
			brother.write_self(file_stream)

			brother = brother.brother
		#endif
	#enddef

	def write_self(self, file_stream):
		if self.child != None:
			file_stream.write('_' + str(self.value) + ' ')
			self.child.write(file_stream)
			file_stream.write(str(self.value) + '_ ')
		else:
			file_stream.write(str(self.value)+ ' ')
		#endif
	#enddef

	def __copy__(self):
		node = TreeNode(copy.copy(self.value))
		node.child = copy.copy(self.child)
		node.brother = copy.copy(self.brother)

		if node.brother != None:
			node.brother.last = node
		#endif

		return node
	#enddef

	@staticmethod
	def from_string(string):
		stripped = string.strip()

		pieces = stripped.split(None, 1)

		if len(pieces) == 2:
			current_string = pieces[1]
		else:
			current_string = ''

		if len(pieces) > 0:
			if len(pieces[0]) > 1 and pieces[0].startswith('_'):
				node = TreeNode(pieces[0][1:])

				node.child, current_string = TreeNode.from_string(current_string)
			elif len(pieces[0]) > 1 and pieces[0].endswith('_'):
				return None, current_string
			else:
				node = TreeNode(pieces[0])

			if current_string != '':
				node.brother, current_string = TreeNode.from_string(current_string)

				if node.brother != None:
					node.brother.last = node
#			else:
#				print("Warning: String ended looking for end of internal node.")
		else:
#			print("Warning: String ended looking for next node.")
			node = None

		return node, current_string
