#! /usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from scipy.ndimage.interpolation import shift
import pyfits
import os.path
from os import mkdir
import numpy as np

class Plotter(object):
    def __init__(self,filelist,step,outdir,ext,clobber):
        self.filelist = filelist
        # datalist holds current state of arrays
        self.datalist = map(pyfits.getdata,filelist)
        self.step = step
        self.outdir = outdir
        self.ext = ext
        self.clobber = clobber
        self.current = 0
        self.orig_data = pyfits.getdata(filelist[0])
        self.active_data = self.datalist[0]
        # store total offsets
        self.offsets = np.zeros((len(filelist),2))
        

        # set up subparser
        self.subparser=argparse.ArgumentParser(description='Parse window text.',prog='')
        self.subparser.add_argument('-s',type=float,help='Step size (default=%.1f)' % self.step)
        self.subparser.add_argument('-w',action='store_true',help="Write current frame to output directory (outdir=%s)"%self.outdir)
        self.subparser.add_argument('-wa',action='store_true',help="Write all frames to output directory (outdir=%s)"%self.outdir)
        self.subparser.add_argument('-wq',action='store_true',help="Write all frames to output directory and quit (outdir=%s)"%self.outdir)
        self.subparser.add_argument('-r',action='store_true',help='Restore original')
        self.subparser.add_argument('--c',action='store_true',help='Force clobber status to True on write')

        self.subparser.add_argument('--q',action='store_true',help='Quit')
        self.fig = plt.figure()
        self.fig.canvas.mpl_disconnect(self.fig.canvas.manager.key_press_handler_id)
        self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.onkey)
        self.pausetext = '-'
        self.pid = None


        # Make directory if doesn't exist
        try:
            mkdir(self.outdir)
        except OSError:
            #directory exists. no big deal.
            pass
        
        print "Type directly into figure"
        print "'left/right/up/down' to translate image"
        print "'</>' to choose previous/next image in input list"
        print "'-h' to see additional options"
        print
        
        #Show initial
        self.display(self.active_data)
        self.displaytext('[0, 0], s=%.2f'%self.step,x=0.60)
        plt.show()
        

    def display(self, data):
        self.fig.clf()
        plt.imshow(data,origin='lower',interpolation='none')
        plt.title(self.filelist[self.current])
        self.fig.canvas.draw()

    def displaytext(self,text,x=0.05,y=0.05,remove=None):
        if remove:
            remove.remove()
        pid = plt.text(x,y,text,color='w',
                       horizontalalignment='left',
                       verticalalignment='bottom',
                       transform=plt.gca().transAxes,
                       path_effects=[PE.withStroke(linewidth=2,foreground='k')])
        self.fig.canvas.draw()
        return pid
        


    def parsetext(self,text):
        args = None
        try:
            # catch -h, or error exit
            args = self.subparser.parse_args(text.split())
        except SystemExit:
            return
            
        if not args:
            return

        if args.c:
            self.clobber = True
            print 'Force clobber on write'
            
        if args.s:
            self.step = args.s
            print 'Step size changed to %.2f' % self.step

        if args.r:
            self.active_data = self.orig_data.copy()
            self.offsets[self.current][0] = 0.0
            self.offsets[self.current][1] = 0.0
            print 'Restored image from %s' % self.filelist[self.current]

        if args.w:
            h = pyfits.getheader(self.filelist[self.current])
            h['N_ORIG_F'] = (self.filelist[self.current],'Original file before nudge')
            h['N_XS'] = (self.offsets[self.current][0],'Xshift of nudge')
            h['N_YS'] = (self.offsets[self.current][1],'Yshift of nudge')
            outfile = os.path.basename(self.filelist[self.current])
            outfile = os.path.splitext(outfile)
            outfile = ''.join([outfile[0],self.ext,outfile[1]])
            outfile = os.path.join(self.outdir,outfile)
            try:
                pyfits.writeto(outfile,data=self.active_data,header=h,clobber=self.clobber)
            except IOError as e:
                print e, "'--c' to force overwrite"
            else:
                print '%s: written to disk, s = [%.2f, %.2f]' % (outfile,self.offsets[self.current][0],self.offsets[self.current][1])

        if args.wq:
            args.wa = True
            args.q = True
                
        if args.wa:
            for idx in range(0,len(self.filelist)):
                h = pyfits.getheader(self.filelist[idx])
                h['N_ORIG_F'] = (self.filelist[idx],'Original file before nudge')
                h['N_XS'] = (self.offsets[idx][0],'Xshift of nudge')
                h['N_YS'] = (self.offsets[idx][1],'Yshift of nudge')
                outfile = os.path.basename(self.filelist[idx])
                outfile = os.path.splitext(outfile)
                outfile = ''.join([outfile[0],self.ext,outfile[1]])
                outfile = os.path.join(self.outdir,outfile)
                try:
                    if idx == self.current:
                        pyfits.writeto(outfile,data=self.active_data,header=h,clobber=self.clobber)
                    else:
                        pyfits.writeto(outfile,data=self.datalist[idx],header=h,clobber=self.clobber)
                except IOError as e:
                    print e, "'--c' to force overwrite"
                    args.q = False #don't quit if file fails to write
                else:
                    print '%s: written to disk, s = [%.2f, %.2f]' % (outfile,self.offsets[idx][0],self.offsets[idx][1])

        if args.q:
            plt.close()
            exit()


        
    def pausekey(self,event):
        if event.key == 'enter':
            self.fig.canvas.mpl_disconnect(self.keycid)
            self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.onkey)
            self.parsetext(self.pausetext)
            self.pausetext = '-'
            self.display(self.active_data)
            self.displaytext('[%.2f, %.2f] s=%.2f'%
                             (self.offsets[self.current][0],
                              self.offsets[self.current][1],
                              self.step),
                             x=0.60)
            return

        elif event.key == 'backspace':
            self.pausetext = self.pausetext[0:-1]
            
        elif len(event.key) > 1:
            return
            
        else:
            self.pausetext = ''.join([self.pausetext,event.key])

        self.pid = self.displaytext(self.pausetext,remove=self.pid)
        return
        
    
    def onkey(self, event):
        if event.key in ['.','>']:
            if self.current >= len(self.filelist)-1:
                return
            self.datalist[self.current] = self.active_data
            self.current += 1
            self.orig_data = pyfits.getdata(self.filelist[self.current])
            self.active_data = self.datalist[self.current]
            
        elif event.key in [',','<']:
            if self.current == 0:
                return
            self.datalist[self.current] = self.active_data
            self.current -= 1
            
            self.orig_data = pyfits.getdata(self.filelist[self.current])
            self.active_data = self.datalist[self.current]

        elif event.key == '-':
             self.fig.canvas.mpl_disconnect(self.keycid)
             self.pausetext = '-'
             self.pid = self.displaytext(self.pausetext)
             self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.pausekey)
             return

        elif event.key == 'left':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[0,-self.step])
            self.offsets[self.current][0] -= self.step

        elif event.key == 'right':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[0,self.step])
            self.offsets[self.current][0] += self.step

        elif event.key == 'down':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[-self.step,0])
            self.offsets[self.current][1] -= self.step

        elif event.key == 'up':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[self.step,0])
            self.offsets[self.current][1] += self.step

        self.display(self.active_data)
        self.displaytext('[%.2f, %.2f] s=%.2f'%
                         (self.offsets[self.current][0],
                          self.offsets[self.current][1],
                          self.step),
                         x=0.60)
        

def main():
    parser = argparse.ArgumentParser(description='Nudge an image, or series of images, by a given step size')
    parser.add_argument('filelist',nargs='+',help='List of input FITS files to be nudged.')
    parser.add_argument('-s',type=float,default=1.0,help='Specify initial step size (default=1.0).')
    parser.add_argument('-o',type=str,default='./nudged/',help="Output directory (default='./nudged/').")
    parser.add_argument('-e',type=str,default='',help="Prefix for filename extension on output (default='').")
    parser.add_argument('--c',action='store_true',help="Clobber on output (default=False.")

    args = parser.parse_args()

    plotter = Plotter(args.filelist,args.s,args.o,args.e,args.c)


if __name__ == '__main__':
    main()







