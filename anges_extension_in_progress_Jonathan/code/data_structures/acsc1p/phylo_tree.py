from decimal import *

# phylogenetic tree (NHX format)
class PhyloTree:
	# PhyloTree constructor
	# root - PhyloTreeNode: the root of the tree
	def __init__(self, root):
		self.root = root

	# finds the ancestor node in the tree
	# Return: PhyloTreeNode - the ancestor node or None if no ancestor node was found
	def find_ancestor(self):
		return self.find_ancestor_rec(self.root)

	# finds the ancestor node below node
	# node - PhylotreeNode: the node to look under
	# Return - PhyloTreeNode: the ancestor node or None if no ancestor node was found
	def find_ancestor_rec(self, node):
		if node.is_ancestor():
			return node

		for child in node.children:
			ancestor = self.find_ancestor_rec(child)

			if ancestor != None:
				return ancestor

		return None

	# returns the number of nodes in the tree
	# Return - int: the number of nodes
	def __len__(self):
		if self.root == None:
			return 0
		else:
			return self.root.num_descendants() + 1

	# returns a string representing the tree in NHX format
	# Return: str - string representation of the tree
	def __str__(self):
		return self.root.subtree_str() + ";"

	# writes the tree to a file in HNX format
	# file_stream: file - file to write to
	def to_file(self, file_stream):
		self.root.to_file(file_stream)
		file_stream.write(';')

	# writes the tree to a file in HNX format
	# file_name: str - name of file to write to
	def to_file_name(self, file_name):
		f = open(file_name, 'w')
		self.to_file(f)
		f.close()

	# Static method which parses a phylogenetic tree from an NHX string
	# string: str - string to parse
	# Return: PhyloTree - the parsed phylogenetic tree (Return.root == None on failure)
	@staticmethod
	def from_string(string):
		node, index = PhyloTreeNode.from_string(string.replace('\r', '').replace('\n', ''))

		return PhyloTree(node)

	# Static method which parses a phylogenetic tree from an NHX file
	# file_name: str - name of the file to read from
	# Return: PhyloTree - the parsed phylogenetic tree (Return.root == None on failure)
	@staticmethod
	def from_file_name(file_name):
		f = open(file_name)
		string = f.read()
		f.close()

		return PhyloTree.from_string(string)

# phylogenic tree node
class PhyloTreeNode:
	# PhyloTreeNode constructor
	# name: str - species name
	# parent: PhyloTreeNode - parent node
	# children: list of PhyloTreeNode - children nodes
	# branch_length: Decimal - length of the branch going from the node to the parent
	# comments: list of str - NHX comments
	def __init__(self, name, parent, children, branch_length, comments):
		self.name = name
		self.parent = parent
		self.children = children
		self.branch_length = branch_length
		self.comments = comments

	# computes the numbe of descendants of the node
	# Return: int - the number of descendants
	def num_descendants(self):
		num = len(self.children)

		for child in self.children:
			num += child.num_descendants()

		return num

	# extracts the species name of the node
	# Return: str - the species name
	def species_name(self):
		for comment in self.comments:
			comment = comment.strip()

			# hard code fix to find ancestor node if it has no species name
			if comment == '@':
				return '@'

			if len(comment) >= 2 and comment.startswith("S=") or comment.startswith("s="):
				return comment[2:].strip()

		return self.name

	# returns True if '@' is contained in the
	# Return - boolean: True if it is the ancestor node and False otherwise
	def is_ancestor(self):
		for comment in self.comments:
			if comment.strip() == '@':
				return True

		return False

	# returns a string representing the tree node in NHX format
	# Return: str - string representation of the tree node
	def __str__(self):
		if len(self.comments) > 0:
			comments_str = '['

			for comment in self.comments:
				comments_str += comment + ':'

			comments_str = comments_str[:-1] + ']'
		else:
			comments_str = ''

		if self.branch_length != Decimal(0):
			branch_str = ":" + str(self.branch_length)
		else:
			branch_str = ""

		return self.name + branch_str + comments_str

	# returns a string representing the sub-tree rooted at this node in NHX format
	# Return: str - string representation of the sub-tree
	def subtree_str(self):
		if len(self.children) > 0:
			string = '('

			for child in self.children:
				string += child.subtree_str() + ','

			string = string[:-1] + ')'
		else:
			string = ''

		return string + str(self)


	# writes the tree node and its children to a file in HNX format
	# file_stream: file - file to write to
	def to_file(self, file_stream):
		if len(self.children) > 0:
			file_stream.write('(')

	 		self.children[0].to_file(file_stream)

			for i in xrange(1, len(self.children)):
				self.children[i].to_file(file_stream)

			file_stream.write(')')

		file_stream.write(str(self))
		file_stream.flush()

	# Method which parses the details (name, branch length, comments) of a phylogenetic tree node from an NHX string
	# string: str - string to parse
	# index: int - index to start parsing in string
	# Return: int - the index that the string was read to (-1 on failure)
	def details_from_string(self, string, index = 0):
		if index >= len(string):
			print("Error: index larger than string length in PhyloTreeNode.details_from_string. string = \"" + string + "\" + index = \"" + str(index) + "\"")

			return -1

		# parse the name of the node
		name_end = PhyloTreeNode.find_first(string, ':[,);', index)
		if name_end == -1:
			print("Warning: cannot find the end of the PhylotreeNode name. string = \"" + string[index:] + "\"")

			return -1

		self.name = string[index:name_end].strip()

		new_index = name_end

		# parse branch length
		if string[new_index] == ':':
			if new_index +1 >= len(string):
				print("Warning: reached end of string looking for branch length. string = \"" + string + "\" index = \"" + new_index + "\"")
				return -1

			branch_end = PhyloTreeNode.find_first(string, "[,);", new_index)
			if branch_end == -1:
				print("Warning: cannot find the end of the PhylotreeNode branch length. string = \"" + string[index:] + "\"")
				return -1

			try:
				self.branch_length = Decimal(string[new_index+1:branch_end].strip())
			except:
				print("Warning: branch_length not a Decimal. branch_length \"" + string[new_index+1:branch_end].strip() + "\"")
				return -1

			new_index = branch_end
		else:
			self.branch_length = Decimal(0)

		# parse comments
		if string[new_index] == '[':
			if new_index + 1 >= len(string):
				print("Warning: reached end of string looking for comments. string = \"" + string + "\" index = \"" + new_index + "\"")
				return -1

			comments_end = string.find(']', new_index)

			if comments_end == -1:
				print("Warning: cannot find the end of the PhylotreeNode comments. string = \"" + string[index:] + "\"")

				return -1

			comments_string = string[new_index + 1:comments_end]

			self.comments = comments_string.split(':')

			new_index = comments_end + 1
		else:
			self.comments = []

		return new_index

	# finds the index of the first occurence of a character in chars in string starting at index
	# string: str - the string to find the characters
	# chars: str - the list of characters to look for
	# index: int - the index in string to start looking at (default = 0)
	# Return: int - the index or -1 if no occurence is found
	@staticmethod
	def find_first(string, chars, index = 0):
		new_index = index

		while new_index < len(string) and not string[new_index] in chars:
			new_index += 1

		if new_index > new_index < len(string):
			return -1

		return new_index

	# Static method which parses a phylogenetic tree node and its children from an NHX string
	# string: str - string to parse
	# index: int - index to start parsing in string
	# parent: PhyloTreeNode - the parent node the node to be parsed
	# Return: PhyloTreeNode, int - the parsed phylogenetic tree node (None on failure), the index of the string that was read to
	@staticmethod
	def from_string(string, index = 0, parent = None):
		new_index = index

		if new_index >= len(string):
			print("Error: index larger than string length in PhyloTreeNode.from_file. string = \"" + string + "\" + index = \"" + str(index) + "\"")
			return None, 0

		node = PhyloTreeNode("", parent, [], Decimal(0), [])

		if string[new_index] == '(':
			while new_index < len(string) and string[new_index] != ')':
				child, new_index = PhyloTreeNode.from_string(string, new_index + 1, node)

				if child == None:
					return None, 0

				node.children.append(child)

			new_index = new_index + 1

			if new_index >= len(string):
				print("Warning: reached end of string looking for ')' in PhyloTreeNode.from_file. string = \"" + string + "\" + index = \"" + str(index) + "\"")
				return None, 0

		new_index = node.details_from_string(string, new_index)

		if new_index == -1:
			return None, 0

		return node, new_index
