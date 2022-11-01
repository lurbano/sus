'''1d horizontal erosion model (diffusional) with vertical columns of soil susceptibility'''
'''  created by Lensyl Urbano 10/16/2009'''
'''  Go to line 250 (or thereabouts) for model parameters'''
'''  Cheers'''
'''   '''
'''  NOTE: Click the left window (titled "VPython") to run the model'''
'''   '''

from visual import *
from visual.graph import *
        

class profile:
    def __init__(self, z, dx, scale=1.0):
        # z is a 1-d array
        # dx is a float
        self.z = z
        self.dx = dx
        self.scale = scale
        self.line = curve(x=self.dx*arange(len(self.z)), y = self.z * self.scale, z=0.0*self.z)
        #print "profile.x=", self.line.x, self.dx

        self.nodes = []
        for i in range(len(self.z)):
            self.nodes.append(sphere(pos=self.line.pos[i], radius=dx/5.0))

    def update(self):
        self.line.y = self.z * self.scale
        for i in range(len(self.z)):
            self.nodes[i].pos=self.line.pos[i]


def mouse_pause(scene_x):
    l_stop = 0
    while l_stop == 0:
        m1 = scene_x.mouse.getevent()
        if m1.click:
            #print "bingo"
            x=1
        if m1.release:
            print "go"
            l_stop = 1

class axes:
    def __init__(self, xmin = 0.0, xmax = 10.0, dx = -1, xscale = 1.0, xcolor = color.red, xshow = True,
                                ymin = 0.0, ymax = 10.0, dy = -1, yscale = 1.0, ycolor = color.blue, yshow = True,
                                zmin = 0.0, zmax = 10.0, dz = -1, zscale = 1.0, zcolor = color.green, zshow = True):
        
        (self.xmin, self.xmax, self.dx, self.xscale, self.xcolor, self.xshow) = (xmin, xmax, dx, xscale, xcolor, xshow)
        (self.ymin, self.ymax, self.dy, self.yscale, self.ycolor, self.yshow) = (ymin, ymax, dy, yscale, ycolor, yshow)
        (self.zmin, self.zmax, self.dz, self.zscale, self.zcolor, self.zshow) = (zmin, zmax, dz, zscale, zcolor, zshow)

        self.nxtics = 11
        self.nytics = 11
        self.nztics = 11
        
        if self.xshow:
            self.draw_xaxis()
        if self.yshow:
            self.draw_yaxis()
        if self.zshow:
            self.draw_zaxis()

    def draw_xaxis(self):
            self.xaxis = curve(pos=[(self.xmin*self.xscale,0), (self.xmax*self.xscale,0)], color = self.xcolor)
            print "xaxis.pos", self.xaxis.pos
            if self.dx == -1:
                self.dx = self.reset_tic_scale(self.xmax, self.xmin, self.nxtics)
            self.tics = []
            for i in range(self.nxtics):
                xpos = (self.xmin+i*self.dx)*self.xscale
                #print xpos
                self.tics.append(tick_mark(pos=vector(xpos, 0, 0), text='%1.0f' % xpos))
            
    def draw_yaxis(self):
            self.zaxis = curve(pos=[(0,self.ymin), (0,self.ymax)], color = self.ycolor)
            if self.dy == -1:
                self.dy = self.reset_tic_scale(self.ymax, self.ymin, self.nytics)
            self.tics = []
            for i in range(self.nytics):
                ypos = self.ymin+i*self.dy
                self.tics.append(tick_mark(pos=vector(0, ypos, 0), text='%1.0f' % ypos, axis="y"))

    def draw_zaxis(self):
            self.zaxis = curve(pos=[(0,0,self.zmin), (0,0,self.zmax)], color = self.zcolor)

    def reset_tic_scale(self, vmax, vmin, ntics):
        return (vmax - vmin) / float(ntics-1)

class tick_mark:
    def __init__(self, pos, text, visible = True, length = 1.0, axis = "x", offset=None, show_tic=True, show_label=True):
        (self.pos, self.text, self.visible, self.length, self.axis) = (pos, text, visible, length, axis)
        if offset == None:
            self.offset = -self.length * 2.0
        if axis == "x":
            self.tic = curve(pos=[(self.pos),(self.pos+vector(0, -self.length, 0))])
            self.label = label(pos=self.pos+vector(0,self.offset,0), text = self.text, box=False, opacity=0)
        if axis == "y":
            self.tic = curve(pos=[(self.pos),(self.pos+vector(-self.length, 0, 0))])
            self.label = label(pos=self.pos+vector(self.offset,0,0), text = self.text, box=False, opacity=0, xoffset=-1, line=0)
            
class sus_plot:
    def __init__(self, z, sus, old_scene, dz=1, show_nodes=False, node_number=-9999):
        self.show_nodes = show_nodes
        title = "Susceptibility Profile" 
        self.graph = display(title=title, x=500, y=0, width=500, height=500)
        self.axes = axes(xmin=0, xmax=1e-6, xscale = 1e7, ymin=-10.0, ymax=0, zshow=False) #xmax=1e-7, xscale = 1e8,
        self.title = label(pos=(0, 2), text="Node: %(n)3i" % {'n' : node_number})
        self.graph.center = vector(self.axes.xscale*(self.axes.xmax-self.axes.xmin)/2.0, \
                                                   self.axes.yscale*(self.axes.ymin-self.axes.ymax)/2.0)
        self.line = curve(x=sus*self.axes.xscale, y=z)
        if self.show_nodes:
            self.nodes = []
            for i in range(len(z)):
                self.nodes.append(sphere(pos=(sus[i]*self.axes.xscale,z[i]), radius=0.1))
        old_scene.select()
    def update(self, z, sus, node_number=-9999):
        #self.graph.title = "Susceptibility: Node %(n)3i" % {'n' : node_number}
        if node_number <> -9999:
            self.title.text="Node: %(n)3i" % {'n' : node_number}
        self.line.x = sus*self.axes.xscale
        self.line.y = z
        if self.show_nodes:
            for i in range(len(self.nodes)):
                self.nodes.pos= vector(sus[i]*self.axes.xscale,z[i])
        #print "line=", max(sus), min(sus)

def update_soil_graph(soilgrf, soilcols, n, slope):
    slope.nodes[soilgrf.node_number].color = color.white
    soilgrf.z = soilcols[n].z
    soilgrf.sus = soilcols[n].sus
    soilgrf.node_number = n
    soilgrf.update_graph()
    slope.nodes[soilgrf.node_number].color = color.red
    
def change_graph_node(s, soilgrf, soilcols, slope):
        if s == "up" or s == "right":
            #print "alpha"
            if soilgrf.node_number < ncells-2:
                n = soilgrf.node_number+1
                update_soil_graph(soilgrf, soilcols, n, slope)
        elif s == "down" or s == "left":
            #print "beta"
            if soilgrf.node_number >= 1:
                n = soilgrf.node_number-1
                update_soil_graph(soilgrf, soilcols, n, slope)

class soil_column:
    def __init__(self, ncells=1000, dz=0.01, node_number=-9999):
        (self.ncells, self.dz, self.node_number) = (ncells, dz, node_number)
        try: #numpy
            self.sus = zeros((ncells,), float)
        except:
            self.sus = zeros((ncells,), Float)
        self.z = arange(0, -ncells*dz, -dz)
        self.dz_top = self.z[0] - self.z[1]

    def add_soil_susceptibility(self, depth, sus_in):
        '''add soil with the given susceptibility (sus_in) to the soil column'''
        ''' for the top cell, do a linear average of the new and old susceptibility'''
        depth = abs(depth)
        #print "depth", depth, self.dz, self.dz_top, self.dz - self.dz_top
        if depth <= (self.dz - self.dz_top): #then add sediment to the upper most cell
            '''calculate the average susceptibility currently in the uppermost cell'''
            sus_old = self.sus[0] #(self.sus[0] + self.sus[1]) / 2.0
            sus_new = (sus_old * self.dz_top + depth * sus_in) / (self.dz_top + depth)
            #print "sus_old, sus_new", sus_old, sus_new
            '''adjust susceptibility of top node given the new average susceptibility (sus_new)'''
            self.sus[0] = sus_new  #2.0 * sus_new - self.sus[1]
            '''adjust location of nodes'''
            self.dz_top += depth
            self.z[1:] -= depth
        else:
            '''fill in the old top node'''
            dz_fill = self.dz - self.dz_top
            sus_old = self.sus[0] # (self.sus[0] + self.sus[1]) / 2.0
            sus_new = (sus_old * self.dz_top + dz_fill * sus_in) / self.dz
            self.sus[0] = sus_new
            self.dz_top = self.dz
            d2 = depth - dz_fill #d2 is the depth remaining after filling top node
            #print "dz_fill, sus_old, sus_new", dz_fill, sus_old, sus_new
            '''add new nodes'''
            nnew_nodes = int(floor(d2 / self.dz) + 1)
            dz_top_new = d2 - self.dz * (nnew_nodes-1)
            self.sus[nnew_nodes:] = self.sus[:self.ncells-nnew_nodes]
            self.sus[:nnew_nodes] = sus_in
            '''recalculate elevations'''
            self.adjust_z(dz_top_new)
        '''check for negative susceptibility'''
        if self.sus[0] < 0.0:
            if depth <= (self.dz - self.dz_top):
                print ""
                x = hendrix #less
            else:
                print ""
                x = hendrix #more
            self.sus[0] = 0.0
                


    def remove_soil_susceptibility(self, depth):
        '''average the susceptibility from the surface (z=0.0) to depth.'''
        ''' use a linear weighted average between nodes'''
        ''' returns the average susceptibility of the removed soil'''
        depth = abs(depth)
        if depth <= self.dz_top:
            zd = 0.5 * depth
            sus_avg = self.sus[0] #(self.sus[0]*(self.dz_top - zd) + self.sus[1] * (zd)) / self.dz_top
            #'''adjust susceptiblity of top node'''
            #self.sus[0] = 2.0 * sus_avg - self.sus[0]
            '''adjust the location of the nodes'''
            self.dz_top = self.dz_top - depth
            self.z[1:] +=  depth
            #print 'self.dz_top', self.dz_top
        else:
            '''use a weighted average with a linear approximation of susceptibility between nodes'''
            ''' first node'''
            sidz = self.dz_top * self.sus[0] # self.dz_top * (self.sus[0] + self.sus[1])/ 2.0
            #print "nodes 0-1", sidz
            n = 2
            ''' all full nodes'''
            while depth > abs(self.z[n]):
                sidz += self.sus[n-1] * self.dz #sidz + ((self.sus[n-1] + self.sus[n]) / 2.0) * self.dz
                #print "nodes", n, sidz
                n += 1
            ''' last node'''
            dz_out = depth + self.z[n-1] #the amount removed from the last node
            #zd = 0.5 * dz_out
            #sus_av = (self.sus[n-1] * (self.dz - zd) + self.sus[n] * zd ) / self.dz
            #sus_newtop = (self.sus[n-1] * (self.dz - dz_out)
            #print "last node", n, self.sus[n-1], self.sus[n], self.dz, zd, sus_av
            sidz += self.sus[n-1] * dz_out #sidz + sus_av * dz_out
            ''' outgoing susceptibility'''
            sus_avg = sidz / depth
            '''adjust susceptibility of top node'''
            self.sus[0] = self.sus[n-1] #2.0 * sus_av - self.sus[n-1]
            '''adjust  the location of the nodes'''
            self.adjust_z(self.dz-dz_out)
            '''adjust the susceptibility of the nodes'''
            self.sus[1:self.ncells-n+1] = self.sus[n:]
            self.sus[self.ncells-n+1:] = 0.0
            #print "indices =", self.ncells-n, n

        if self.dz_top == 0.0:
            '''shift nodes downward'''
            #print "shift"
            self.z[1:-1] = self.z[2:]
            self.z[-1] -= self.dz
            self.dz_top = 0.0
            self.sus[1:-1] = self.sus[2:]
            self.sus[-1] = 0.0

        return sus_avg

    def adjust_z(self, dz_top_new):
        self.dz_top = dz_top_new
        self.z[1:] = arange(0.0, -(self.ncells-1)*self.dz, -self.dz) - self.dz_top
            
    def create_graph(self):
        self.graph = sus_plot(self.z, self.sus, scene, show_nodes=False, node_number = self.node_number)

    def update_graph(self):
        self.graph.update(self.z, self.sus, node_number=self.node_number)

def read_topo_from_output(datafile="suscept.txt"):
    inf = open(datafile, 'r')
    ln = inf.readline()
    ln = inf.readline()
    elev = ln.split(",")
    elev = elev[2:-1]
    for i in range(len(elev)):
        elev[i] = float(elev[i])
    print elev[0], len(elev)
    return array(elev)    
    
'''***************************************'''
'''time dimensions'''
'''***************************************'''
dt = 10.0
nsteps = 1001
print_out_step = 10 #print every so many steps


'''***************************************'''
'''Erosion Model parameters'''
'''***************************************'''
k = 0.01         #diffusion coefficient
rho = 0.8       #bulk density of the soil


'''***************************************'''
'''Model physical dimensions'''
'''***************************************'''
ncells = 100
dx = 1.0
'''initial topography'''
# set up initial arrays
try: #numpy
    topo = ones((ncells,), float)*50.0
except:
    topo = ones((ncells,), Float)*50.0
topo[ncells/2:ncells] = 25.0
# input initial topgraphy (two options)
#  option 1: copy and paste values
#topo = array([49.9521499317,49.9519271202,49.9496931132,49.9453889506,49.938902203,49.9300662329,49.9186591538,49.9044025494,49.8869600247,49.8659356835,49.8408726397,49.8112516845,49.7764902438,49.7359417715,49.6888957283,49.6345783021,49.5721540256,49.5007284386,49.4193519379,49.3270249384,49.2227044513,49.1053121585,48.9737440306,48.8268815005,48.6636041632,48.4828039266,48.2834004929,48.064357999,47.8247025964,47.5635407037,47.2800776193,46.9736361449,46.6436748351,46.289805467,45.9118093077,45.5096517554,45.0834949403,44.6337078897,44.1608739013,43.6657948109,43.1494919049,42.6132032907,42.0583776236,41.4866641655,40.8998992462,40.3000892819,39.6893905982,39.070086387,38.4445612022,37.8152734676,37.1847265251,36.5554387894,35.9299136022,35.3106093871,34.6999106974,34.100100724,33.5133357917,32.941622315,32.3867966211,31.8505079691,31.3342050095,30.8391258439,30.3662917499,29.9165045522,29.4903475327,29.0881896976,28.7101931485,28.3563232455,28.0263612047,27.7199187353,27.436454302,27.1752905882,26.9356327372,26.7165869651,26.51717916,26.3363731183,26.173088104,26.0262154633,25.8946340746,25.7772244618,25.6728814465,25.5805252667,25.4991111266,25.4276371918,25.365151072,25.3107548717,25.2636089095,25.2229342308,25.1880140526,25.158194287,25.1328832984,25.1115510475,25.0937277735,25.0790023586,25.0670205111,25.0574828909,25.05014329,25.0448069638,25.0413291968,25.0396141684])
#  option 2: read from existing output file (default "suscept.txt" but this file will be overwritten at the end of this simulation)
#topo = read_topo_from_output("suscept.txt")
print "topo", topo

'''***************************************'''
'''Set up soil columns for magnetic susceptibility'''
'''***************************************'''
n_soil_cells = 1000 # for 10 m at dz = 0.01
soil_dz = 0.01 # 1cm
chi_surface = 1.0e-9
z_star = 1.4

'''initial soil depths at nodes'''
z = arange(0, -n_soil_cells*soil_dz,-soil_dz) #depth below the surface
soilcols = []
for i in range(ncells):
    soilcols.append(soil_column(n_soil_cells, soil_dz, node_number=i))

'''***************************************'''
'''Output'''
'''***************************************'''
'''output column number'''
outcol = 65
'''output file name'''
outfile = "suscept.txt"

'''create output graphs'''
soilgrf = soilcols[outcol]
soilgrf.create_graph()

'''center output window'''
try: #numpy
    scene.center = vector(ncells*dx/2.0,(amax(topo)+amin(topo))/2.0,0)
except:
    scene.center = vector(ncells*dx/2.0,(max(topo)+min(topo))/2.0,0)
    
'''slope topography'''
slope = profile(topo, dx)
slope.nodes[outcol].color = color.red


'''***************************************'''
'''model setup'''
'''***************************************'''
'''storage arrays'''
try: #numpy
    q_s = zeros(topo.shape, float)
    dh = zeros(topo.shape, float)
    dh_in = zeros(topo.shape, float)
    dh_out = zeros(topo.shape, float)
    topo_old = zeros(topo.shape, float)
    sus_out = zeros(topo.shape, float)
except:
    q_s = zeros(topo.shape, Float)
    dh = zeros(topo.shape, Float)
    dh_in = zeros(topo.shape, Float)
    dh_out = zeros(topo.shape, Float)
    topo_old = zeros(topo.shape, Float)
    sus_out = zeros(topo.shape, Float)

'''labels'''
try: #numpy
    timestamp = label(pos=((ncells*dx*0.8, amax(topo)*1.2, 0)))
    timestamp.text = "Timestep = " + str(0)
    xyaxes = axes(xmax=ncells*dx, ymax=amax(topo)*1.1, zshow=False)
except: #Numeric
    timestamp = label(pos=((ncells*dx*0.8, max(topo)*1.2, 0)))
    timestamp.text = "Timestep = " + str(0)
    xyaxes = axes(xmax=ncells*dx, ymax=max(topo)*1.1, zshow=False)

    
print "\n Ready to run simulation. \n Click on VPython window to start simulation \n"
mouse_pause(scene)

'''***************************************'''
'''Timesteps'''
'''***************************************'''
print "Running simulation"
for i in range(nsteps):

    ''' q_s = forward sediment flux (in the x direction) (flow out)'''
    '''diffusion'''
    q_s[:-1] = k * (topo[1:] - topo[:-1])/dx

    '''dh = change in topography for internal cells'''
    dh_in[1:] = -q_s[:-1] * dt / rho
    dh_out[:-1] =  dh_in[1:]  
    dh = dh_in - dh_out

    '''store old topo in case of emergency'''
    topo_old[:] = topo[:]

    '''boundary conditions'''
    '''  right boundary'''
    # bc = constant elevation: then simply don't change the elevation at the boundaries and there will be no sediment flux (q_s)
    #  so do nothing for topo[ncells]
    # NO-SEDIMENT FLUX THEN USE
    #topo[-2] = topo[-1]
    
    '''  left boundary'''
    # bc = no flow on left (x=0): then make the slope on the left 0.0
    topo[0] = topo[1]

    '''adjust topography'''
    topo = topo + dh

    '''update soil profiles'''
    ''' upstream boundary (leftmost node) ages in the usual way'''
    soilcols[0].sus += chi_surface * exp(soilcols[0].z/z_star)
    for n in range(len(soilcols)-1):
        '''advective transport of sediment's of pre-existing susceptibility
            as with the rest of the model, we use a backward looking (timewise) explicit function'''
        '''remove eroded sediment'''
        sus_out[n] = soilcols[n].remove_soil_susceptibility(dh_out[n])
        '''deposit sediment/susceptibility'''
        soilcols[n].add_soil_susceptibility(dh_in[n], sus_out[n-1])
        '''susceptibility pedogenesis function'''
        soilcols[n].sus += chi_surface * exp(soilcols[n].z/z_star)

    '''ouput during simulation'''
    if i % print_out_step == 0:
        print "Timestep=", i
        print "  Time = ", i*dt
        timestamp.text = "Timestep = " + str(i)
        slope.z = topo
        slope.update()
##        print 'dh, topo_old[outcol-1], topo_old[outcol]', dh[outcol]
##        print  topo_old[outcol-1], topo_old[outcol], topo_old[outcol+1]
##        print "sus_in, sus_out", sus_out[outcol-1], sus_out[outcol]
##        print "sus", soilcols[outcol-1].sus[0],  soilcols[outcol].sus[0],  soilcols[outcol+1].sus[0]
        #soilscene.update(z, sus[40])
        #soilcols[outcol].update_graph()
        soilgrf.update_graph()


'''***************************************'''
'''Write output'''
'''***************************************'''
print "Writing output to file:", outfile
outf = open(outfile, "w")
'''header'''
outln = "Node, Depth, "
for i in range(len(soilcols)):
    outln += str(i)+","
outln += "\n"
outf.write(outln)

'''write surface elevations of nodes'''
outln = "Elevation, ,"
for n in slope.z:
    outln += str(n) + ","
outln += "\n"
outf.write(outln)


'''main output table'''
for r in range(len(soilcols[0].sus)):
    outln = str(r) + "," + str(r*soil_dz) + ","
    for c in soilcols:
        outln += "%1.5e," % c.sus[r]
    outln += " \n"
    outf.write(outln)
outf.close()
print "Done writing"

'''***************************************'''
'''Post simulation data analysis'''
'''***************************************'''
print ""
print "Use arrow keys (or click on node) to select nodes to view on graph. \n"
while 1:
    if scene.mouse.events:
        m1 = scene.mouse.getevent()
        if m1.click and m1.release:
            for n in range(len(slope.nodes)):
                if m1.pick == slope.nodes[n]:
                    print "Node selected =", n
                    update_soil_graph(soilgrf, soilcols, n, slope)

    
    if soilgrf.graph.graph.kb.keys: # is there an event waiting to be processed?
        s = soilgrf.graph.graph.kb.getkey() # obtain keyboard information
        #print s
        change_graph_node(s,  soilgrf, soilcols, slope)
            
    if scene.kb.keys: # is there an event waiting to be processed?
        s = scene.kb.getkey() # obtain keyboard information
        #print s
        change_graph_node(s,  soilgrf, soilcols, slope)
            
    

