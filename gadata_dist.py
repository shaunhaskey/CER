import MDSplus 
import numpy
import time
import sys

class gadata_dist:
	"""GA Data Obj"""
	def __init__(self,signal,shot,tree=None,nomds=False):

		# Save object values
		self.signal 		= signal
		self.shot   		= shot
		self.zdata              = -1
		self.xdata              = -1
		self.ydata 		= -1
		self.zunits             = ''
                self.xunits             = ''
		self.yunits		= ''
		self.rank 		= -1

		## Retrieve Data 
                t0 =  time.time()
		found = 0

		# Retrieve data from MDSplus (dist)
  		if nomds == False:
   			try:     
				if tree != None:
					tag 	= self.signal
					fstree 	= tree 
				else:
					# WARNING - This TDI function is time consuming over distributed MDSplus
					#           Performance of Thin:Dist is 20x (holds true for IDL as well)  
					#	    This is due entirely to MDSopen(D3D)
					tag 		= MDSplus.Data.execute('findsig("'+self.signal+'",_fstree)').value
  					fstree    	= MDSplus.Data.execute('_fstree').value
 
				tree = MDSplus.Tree(fstree,shot)
				node = tree.getNode(tag)

				self.zdata  	= node.getData().data()
				self.zunits 	= node.getUnits().value
				self.rank	= numpy.rank(self.zdata)
				self.xdata	= node.getDimensionAt(0).data()
				self.xunits	= node.getDimensionAt(0).getUnits().value
				if self.rank > 1:
					self.ydata      = node.getDimensionAt(1).data()
       		                        self.yunits     = node.getDimensionAt(1).getUnits().value
	
				found = 1	

				# MDSplus seems to return 2-D arrays transposed.  Change them back.
				if numpy.rank(self.zdata) == 2: self.zdata = numpy.transpose(self.zdata)
				if numpy.rank(self.ydata) == 2: self.ydata = numpy.transpose(self.ydata)
				if numpy.rank(self.xdata) == 2: self.xdata = numpy.transpose(self.xdata)

                	except Exception,e:
				print '   Signal not in MDSplus: %s' % (signal,) 
		

		# Retrieve data from PTDATA
                if found == 0:
                        self.zdata = MDSplus.Data.execute('_s = ptdata2("'+signal+'",'+str(shot)+')')
                        if len(self.zdata) != 1:
                                self.xdata = MDSplus.Data.execute('dim_of(_s)')
                                self.rank = 1
                                found = 1

		# Retrieve data from Pseudo-pointname 
		if found == 0:
			self.zdata = MDSplus.Data.execute('_s = pseudo("'+signal+'",'+str(shot)+')')
			if len(self.zdata) != 1:
				self.xdata = MDSplus.Data.execute('dim_of(_s)')
				self.rank = 1
			        found = 1

		if found == 0: 
			print "   No such signal: %s" % (signal,)
			return

                print '   GADATA Retrieval Time : ',time.time() - t0


