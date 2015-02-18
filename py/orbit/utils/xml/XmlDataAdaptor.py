"""
XmlDataAdaptor is a DataAdaptor that specifically supports (reading/writing)
(from/to) XML. XmlDataAdaptor uses xml.dom.minidom package.  In particular, 
you can use the methods of DataAdaptor for all data manipulation.  

To create a new, empty XML document simply invoke:
<code> document_adaptor = XmlDataAdaptor(name); </code>

You can populate the document by adding child nodes.  You can only add a single node to the 
top document node, but otherwise you can add and nest as many nodes as needed. Each 
such node is returned as a DataAdaptor.  For example, to add a node to the top 
document node, invoke:
<code> childAdaptor = document_adaptor.createChild("nodeName") </code>

You can set attributes of nodes with some basic types such as boolean, integer, double
and string where you supply the name of the attribute and the value.  If you add an Object
as a value, then toString() is invoked to fetch the value as a string.  Some examples:
<code> adaptor.setValue("attribute", "some string"); </code>
<code> adaptor.setValue("attr2", 3.14); </code>

You can write an XML document to a file.
For example, to write an XML document to a file invoke:
<code> document_adaptor.writeToFile("file_path"); </code>

You can read an XML document from a string of XML text, a URL or a file.  You may 
specify whether or not to use DTD validation.
For example, to read an XML document from a file invoke:
<code> document_adaptor = XmlDataAdaptor.adaptorForFile("file_path"); </code>

You can fetch named child nodes from a parent node.  Some examples are:
<code> List xAdaptors = parentAdaptor.childAdaptors("X") </code>
<code> List allAdaptors = parentAdaptor.childAdaptors() </code>
<code> DataAdaptor yAdaptor = parentAdaptor.childAdaptor("Y") </code>

You can test if a node defines an attribute:
<code> boolean_status = adaptor.hasAttribute("attribute"); </code>

You can read the value of an attribute:
<code> double_value = adaptor.doubleValue("attr2"); </code>

There are several more methods available for DataAdaptor and XmlDataAdaptor, but 
the examples listed above should provide an overview and foundation.

@author  shishlo
"""

import types

#import python XML DOM parser
import xml.dom.minidom

from orbit.utils import NamedObject, ParamsDictObject

class XmlDataAdaptor(NamedObject,ParamsDictObject):
	def __init__(self, name = "root"):
		"""
		The XmlDataAdaptor constructor. The default name is 'root'. 
		It is an analog of the XAL XmlDataAdaptor class. The Attributes are replaced by parameters.
		"""
		NamedObject.__init__(self,name)
		ParamsDictObject.__init__(self)	
		self.data_adaptors = []
		
	def clean(self):
		""" Remove all data. """
		ParamsDictObject.__init__(self)	
		self.data_adaptors = []	
		
	def hasAttribute(self,attribute):
		""" Returns true or false if the node has the specified attribute """
		return self.hasParam(attribute)
		
	def stringValue(self,attribute):
		""" string value associated with the specified attribute """
		return str(self.getParam(attribute))
		
	def doubleValue(self,attribute):
		""" double value associated with the specified attribute """
		return float(self.getParam(attribute))
		
	def intValue(self,attribute):
		""" integer value associated with the specified attribute """
		return int(self.getParam(attribute))
		
	def booleanValue(self,attribute):
		""" boolean value associated with the specified attribute """
		return bool(self.getParam(attribute))
	
	def doubleArrayValue(self,attribute):
		"""  array of doubles associated with the specified attribute """
		arr = []
		st = self.stringValue(attribute)
		st_arr = st.split()
		for s in st_arr:
			arr.append(float(s))
		return arr
	
	def intArrayValue(self,attribute):
		"""  array of integers associated with the specified attribute """
		arr = []
		st = self.stringValue(attribute)
		st_arr = st.split()
		for s in st_arr:
			arr.append(int(s))
		return arr	
		
	def intArrayValue(self,attribute):
		"""  array of booleans associated with the specified attribute """
		arr = []
		st = self.stringValue(attribute)
		st_arr = st.split()
		for s in st_arr:
			arr.append(bools(s))
		return arr
		
	def setValue(self,attribute,value):
		""" set the value for the attribute """
		if(isinstance(value,(types.TupleType,types.ListType))):
			st = ""
			for s in value:
				st = st+" "+str(s)
			self.setParam(attribute,st)
			return
		self.setParam(attribute,str(value))
		
	def nodeCount(self):
		""" return the number of child node adaptors """
		return len(self.data_adaptors)
		
	def childAdaptors(self,name = None):
		""" return the array of child node adaptors. All of them or with particular name"""
		if(name == None):
			return self.data_adaptors
		adaptors = []
		for adaptor in self.data_adaptors:
			if(adaptor.getName() == name):
				adaptors.append(adaptor)
		return adaptors
		
	def createChild(self,name):
		""" create and return the new child data adaptor """
		adaptor = XmlDataAdaptor(name)
		self.data_adaptors.append(adaptor)
		return adaptor
		
	def attributes(self):
		""" return the list of attributes """
		return self.keys()
		
	def removeAttribute(self,attribute):
		""" remove the particular attribute """
		if(self.hasAttribute(attribute)):
			del self.getParamsDict()[attribute]
		
	def writeToFile(self,file_name):
		""" write the xml document to the file """
		doc = xml.dom.minidom.Document()
		root = doc.createElement("xml_root_node")
		self.makeDomElement(doc,root,self)
		doc.appendChild(root)
		fl_out = open(file_name,"w")
		xml_text = doc.toprettyxml(indent=" ")
		fl_out.write(xml_text)
		fl_out.close()
		
	def makeDomElement(self,doc,parent_element,adaptor):
		""" generate XML tree element to write it to the file """
		element = doc.createElement(adaptor.getName())
		for attribute in adaptor.attributes():
			value = adaptor.stringValue(attribute)
			element.setAttribute(attribute,value)
		for adaptor_child in adaptor.data_adaptors:
			self.makeDomElement(doc,element,adaptor_child)
		parent_element.appendChild(element)
		
		
			

