"""
Author: David Moran, email: drm285@cornell.edu
Supervisor: Professor Hugh Gauch, email: hgg1@cornell.edu
"""
from os import remove, startfile, path, getenv, chdir
from re import search
from tkFileDialog import askopenfilename
import Tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from PIL import Image, ImageTk
from time import sleep

def enable(ent):
    """Enable a Tkinter widget."""
    ent.config(state=tk.NORMAL)
    ent.update()

def disable(ent):
    """Disable a Tkinter widget."""
    ent.config(state=tk.DISABLED)
    ent.update()

def imtolab(image):
    """Return Tkinter photo object from jpg image."""
    image = Image.open(image)
    return ImageTk.PhotoImage(image)

class Terms(object):
    """Terms class represents the terms and conditions window."""

    def __init__(self, master):
        """Launch Tkinter terms and conditions window."""
        self.master = master
        master.protocol('WM_DELETE_WINDOW', self._quit)
        master.title(string="AMMISOFT: Terms and Conditions")
        terms = imtolab('termsconditions.jpg')
        can = tk.Canvas(master, width=450, height=500)
        can.pack(side='top')
        can.create_image(230, 250, image=terms)
        lic = tk.Button(master, text="License", command=lambda:
                        startfile('license.txt'))
        agree = tk.Button(master, text="Agree", command=self.launch_start)
        can.create_window(320, 490, anchor='s', window=lic)
        can.create_window(380, 490, anchor='s', window=agree)
        master.mainloop()

    def _quit(self):
        """Destroy terms and conditions window."""
        self.master.quit()
        self.master.destroy()

    def launch_start(self):
        """Destroy terms window and launch start window."""
        self._quit()
        master = tk.Tk()
        Start(master)
        master.mainloop()

class Start(object):
    """Start class represents the start window."""

    def __init__(self, master):
        """Launch Tkinter start window."""
        self.master = master
        master.title("AMMISOFT: Main Menu")
        master.minsize(width=400, height=350)
        master.maxsize(width=500, height=450)
        master.protocol('WM_DELETE_WINDOW', self._quit)
        self.output = None
        tk.Button(master, text="Choose File", command=self.browse).pack(side='top')
        self.path = tk.Text(master, wrap=tk.WORD, width=40, height=2,
                            borderwidth=4, font=('timesnewroman', 10))
        self.path.pack(side='top')
        self.fort = tk.Button(master, text="AMMI Analysis", command=self.openammi)
        self.fort.pack(side='top')
        disable(self.fort)
        labcan = tk.Canvas(master, width=150, height=110)
        self.dirlab = {}
        self.dirlab['choose'] = tk.Label(labcan, text="Please choose a valid AMMI file.")
        self.dirlab['input'] = tk.Label(labcan, text=("AMMI input file chosen." + 
                                                      "\n" + "Ready to perform analysis."))
        self.dirlab['output'] = tk.Label(labcan, text=("AMMI output file chosen." + 
                                                       "\n" + "Ready for graphing."))
        with open('ammisuccess.txt', 'r') as myfile:
            self.dirlab['success'] = tk.Label(labcan, text=myfile.read())
        self.dirlab['choose'].pack(side='top')
        labcan.pack(side='top')
        self.tab = tk.Button(master, text="View Tables",
                             command=lambda: startfile(self.output))
        self.tab.pack(side='top')
        disable(self.tab)
        self.run = tk.Button(master, text="Make Graphs", command=self.launch_graph)
        self.run.pack(side='top')
        disable(self.run)
        can = tk.Canvas(master, width=150, height=10)
        tk.Button(can, text="Data Format", command=lambda:
                  startfile('AMMISOFTDataFormatNotepad.txt')).pack(side='left')
        tk.Button(can, text="Data Format Short Guide", command=lambda:
                  startfile('AMMISOFTDataFormatShortGuide.txt')).pack(side='left')
        can.pack(side='bottom')

    def _quit(self):
        """Destroy start window."""
        self.master.quit()
        self.master.destroy()

    def launch_graph(self):
        """Destroy start window and launch graphing window"""
        self._quit()
        master = tk.Tk()
        Graph(master, self.output)
        master.mainloop()

    def ammirun(self):
        """Run Fortran program AMMIFORT.exe."""
        sfpath = 'InHouseSuccessFailure.txt'
        if path.exists(sfpath):
            remove(sfpath)
        sleep(.2)
        startfile('AMMIFORT.exe')
        sleep(.5)
        if path.exists(sfpath):
            errorfile = open(sfpath)
            msgs = errorfile.read()
            if msgs[0:7] == "Success":
                disable(self.fort)
                enable(self.tab)
                enable(self.run)
                for lab in self.dirlab.keys():
                    self.dirlab[lab].pack_forget()
                self.dirlab['success'].pack(side='top')
            else:
                root = tk.Tk()
                def okay():
                    root.quit()
                    root.destroy()
                    remove(self.output)
                tk.Label(root, text=msgs).pack(side='top')
                tk.Button(root, text="OK", command=okay).pack(side='top')
                root.mainloop()

    def browse(self):
        """Launch browsing window for input file."""
        stpath = path.join(getenv('USERPROFILE'), "Documents\AMMISOFT")
        filename = askopenfilename(parent=self.master,
                                   title='Choose a valid AMMISOFT file',
                                   filetypes = [('text files', '.txt')],
                                   initialdir = stpath)
        open('InHouseFilenames.txt', 'w').close()
        if filename[-8:] == 'ammi.txt' or filename[-4:] == '.txt':
            for lab in self.dirlab.keys():
                self.dirlab[lab].pack_forget()
            self.path.delete('1.0', tk.END)
            self.path.insert(tk.END, filename)
            self.path.update()
            if filename[-8:] == 'ammi.txt':
                disable(self.fort)
                enable(self.tab)
                enable(self.run)
                self.output = filename
                self.dirlab['output'].pack(side='top')
            elif filename[-4:] == '.txt':
                enable(self.fort)
                disable(self.tab)
                disable(self.run)
                readfile = open('InHouseFilenames.txt', 'w+', 0)
                self.output = filename[:-4] + 'ammi' + filename[-4:]
                readfile.write(filename + '\n' + self.output)
                self.dirlab['input'].pack(side='top')

    def openammi(self):
        """Protocol for when the 'AMMI Analysis' button is clicked."""
        if path.exists(self.output):
            root = tk.Tk()
            def _remove():
                remove(self.output)
                root.quit()     # stops mainloop
                root.destroy()
                self.ammirun()
            def dontremove():
                root.quit()
                root.destroy()
            tk.Label(root, text=('AMMI analysis has already been done ' +
                                 'for this dataset.' + '\n' +
                                 'Overwrite previous analysis?')).pack(side='top')
            yesno = tk.Canvas(root, width=150, height=10)
            tk.Button(yesno, text="Yes", command=_remove).pack(side='left')
            tk.Button(yesno, text="No", command=dontremove).pack(side='right')
            yesno.pack(side=tk.TOP)
            root.mainloop()
        else:
            self.ammirun()

class Graph(object):
    """Graph class represents the graphing window."""

    def __init__(self, master, output):
        """Read data and launch Tkinter graphing window."""
        self.master = master
        master.title("AMMIGRAPH")
        master.state('zoomed')
        master.protocol('WM_DELETE_WINDOW', self.launch_start)
        f = Figure(dpi=110)
        self.a = f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(f, master=master)
        self.canvas.show()
        self.higher = []
        self.gen, self.env, self.mean, self.higher = self.readdata(output)
        self.x_winner, self.y_winner, self.y_min, self.y_max, self.winnerindex = self.linedata()
        self.winners = False
        self.infopic = {}
        modkey = ['ammi1', 'linear', 'megaenv', 'higher']

        for pic in modkey:
            self.infopic[pic] = imtolab(pic + 'story.jpg')
            
        self.pic = self.infopic['ammi1']
        self.color = ['blue', 'red', 'b^', 'ro']
        self.exitflag = True

        self.model = tk.StringVar()
        bcanvas = tk.Canvas(master, width=150, height=900)
        
        for key, label in zip(modkey, ("AMMI1", "Linear model",
                                       "AMMI1 winners", "AMMI2 or higher")):
            r = tk.Radiobutton(bcanvas, text=label, variable=self.model,
                               value=key, command=self.biplot)
            r.pack(side='top', anchor='w')
            if key == 'ammi1':
                r.select()
        tk.Button(bcanvas, text="Show graph explanation",
                  command=self.info).pack(side='top')
        self.createline(bcanvas)
        self.chkbox = []

        for var, marklab in zip(('mark', 'lab'), (' markers', ' labels')):
            for meas, lab in zip((self.gen, self.env),
                                 ('Genotype', 'Environment')):
                meas[var] = tk.IntVar()
                c = tk.Checkbutton(bcanvas, text=lab + marklab,
                                   variable=meas[var], command=self.biplot)
                self.chkbox.append(c)
                c.pack(side=tk.TOP, anchor=tk.NW)
                c.select()
            self.createline(bcanvas)

        self.col = tk.StringVar()
        
        for val, lab in zip(('col', 'bw'), ("Color", "Black/White")):
            r = tk.Radiobutton(bcanvas, text=lab, variable=self.col,
                               value=val, command=self.biplot)
            r.pack(side='top', anchor='nw')
            if val == 'col':
                r.select()
        self.createline(bcanvas)

        self.highcan = tk.Canvas(bcanvas, width=150, height=10)
        labcan = tk.Canvas(self.highcan, width=150, height=10)
        
        for lab in ('x-axis', 'y-axis'):
            tk.Label(labcan, text=lab).pack(side='left')
            
        labcan.pack(side='top')
        self.xvar = tk.StringVar()
        self.yvar = tk.StringVar()
        self.defhigh = []

        for lev in self.higher:
            canvas = tk.Canvas(self.highcan, width=150, height=10)
            for var in (self.xvar, self.yvar):
                r = tk.Radiobutton(canvas, text=lev, variable=var, value=lev,
                                   command=self.biplot)
                if ((lev == "IPC1" and var == self.xvar) or
                        (lev == "IPC2" and var == self.yvar)):
                    r.select()
                r.pack(side='left')
            canvas.pack(side='top')
            
        bcanvas.pack(side='left', anchor='nw')
        toolbar = NavigationToolbar2TkAgg(self.canvas, master)
        toolbar.pack(side='bottom', anchor='e')
        self.biplot()
        self.canvas._tkcanvas.pack(side='top', fill='both', expand=True)

    def launch_start(self):
        """Destroy graphing window and launch start window."""
        self.master.quit()
        self.master.destroy()
        master = tk.Tk()
        Start(master)
        master.mainloop()

    def readdata(self, output):
        """Read data from AMMIFORT output file."""
        gen = {}
        env = {}
        ammifile = open(output)
        
        for line in ammifile:
            line = line.rstrip()
            if search("Computer-readable AMMI Parameters", line):
                pstring = next(ammifile)
                gen['Num'] = int(pstring[1:5])
                env['Num'] = int(pstring[6:10])
                mean = float(pstring[21:40])
                lev = int(pstring[16:20])
                higher = ['IPC' + str(k + 1) for k in range(lev)]
                          
                for meas in (gen, env):
                    for par in ['Name', 'Mean'] + higher:
                        meas[par] = []

                    for x in range(meas['Num']):
                        dstring = next(ammifile)
                        meas['Name'].append(dstring[7:11])
                        meas['Mean'].append(float(dstring[12:31]))
                        beg = 12
                        end = 31
                        for par in higher:
                            beg += 20
                            end += 20
                            meas[par].append(float(dstring[beg:end]))
                            
        return gen, env, mean, higher

    def linedata(self):
        """Calculate variables for lineplot."""
        winnerindex = []
        x_winner = [min(self.env['IPC1']), max(self.env['IPC1'])]
        y_winner = []
        y_min = float('inf')
        
        for minmax in (x_winner[0], x_winner[1]):
            y_max = float('-inf')
            
            for i in range(self.gen['Num']):
                y_test = self.gen['Mean'][i] + self.gen['IPC1'][i] * minmax
                if y_test > y_max:
                    y_max = y_test
                    index = i
                if y_test < y_min:
                    y_min = y_test
                    
            winnerindex.append(index)
            y_winner.append(y_max)
            
        y_max = y_winner[0] if y_winner[0] > y_winner[1] else y_winner[1]
        x_bound = x_winner[0]
        index1 = winnerindex[0]
        check1 = False
        
        while index1 != winnerindex[1]:
            check2 = False
            x = x_winner[1]

            for i in range(0, self.gen['Num']):
                if i != index1 and (self.gen['IPC1'][i] - self.gen['IPC1'][index1]) != 0:
                    x_test = (self.gen['Mean'][index1] -
                              self.gen['Mean'][i])/(self.gen['IPC1'][i] -
                                                    self.gen['IPC1'][index1])
                    if x_test < x and x_test > x_bound:
                        x = x_test
                        y = self.gen['Mean'][i] + self.gen['IPC1'][i] * x
                        index = i
                        check1 = check2 = True
                        
            if check2:
                winnerindex.append(index)
                x_winner.append(x)
                y_winner.append(y)
                x_bound = x
                index1 = index
            else:
                break
        if check1:
            winnerindex.pop()
        return x_winner, y_winner, y_min, y_max, winnerindex

    def createline(self, can):
        canvas1 = tk.Canvas(can, width=150, height=10)
        canvas1.create_line(0, 10, 150, 10)
        canvas1.pack(side='top', anchor='w')

    def info(self):
        """Show information window for selected graph model."""
        self.exitflag = False
        self.infwin = tk.Toplevel()
        self.infwin.title('Graph Explanation')
        self.infwin.protocol('WM_DELETE_WINDOW', self.infoquit)
        label = tk.Label(self.infwin, image=self.pic)
        label.configure(background='white')
        label.pack()
        self.infwin.mainloop()

    def infoquit(self):
        """Destroy information window."""
        self.exitflag = True
        self.infwin.quit()
        self.infwin.destroy()

    def biplot(self):
        """Plot a new biplot, or call the lineplot method."""
        if not self.exitflag:
            self.infoquit()
            
        self.pic = self.infopic[self.model.get()]
        self.color = (['black', 'black', 'k^', 'ko'] if self.col.get() == 'bw'
                      else ['blue', 'red', 'b^', 'ro'])
        model = [m == self.model.get() for m in ('ammi1', 'megaenv', 'higher')]
                 
        if model[0] or model[2]:
            for c in self.chkbox:
                enable(c)
        else:
            for c in self.chkbox:
                disable(c)
                
        if model[0] or model[1]:
            xlabel = 'Mean'
            ylabel = 'IPC1'        
            x1 = self.gen['Mean']
            y1 = self.gen['IPC1']
            x2 = self.env['Mean']
            y2 = self.env['IPC1']

            self.highcan.pack_forget()
        elif model[2]:
            self.highcan.pack(side='top')
            xlabel = self.xvar.get()
            ylabel = self.yvar.get()
            x1 = self.gen[xlabel]
            y1 = self.gen[ylabel]
            x2 = self.env[xlabel]
            y2 = self.env[ylabel]
        else:
            self.lineplot()
            return None
            
        a = self.a
        a.clear()
        a.set_aspect('auto')
        x3 = []
        y3 = []

        for da, lim2, out in zip((x1+x2, y1+y2), (.12000, .07000), (x3, y3)):
            drange = abs(max(da)-min(da))
            for co in (.030000, lim2, .0050000, .0010000):
                out.append(co * drange)
            out[0] = min(da) - out[0]
            out[1] = max(da) + out[1]

        a.set_xlabel(xlabel, fontsize=20)
        a.set_ylabel(ylabel, fontsize=20)
        a.set_xlim(x3[0], x3[1])
        a.set_ylim(y3[0], y3[1])
        fnt = 12
        
        if model[0]:
            a.axhline(y=0, color='k')
            a.axvline(x=self.mean, color='k')
        elif model[2]:
            a.axhline(y=0, color='k')
            a.axvline(x=0, color='k')
            a.set_aspect('equal', adjustable='box')
        else:
            fnt = 17
            for line in zip(self.x_winner[2:]):
                a.axhline(y=line, color='k', linewidth=2)
                
        for meas, xl, yl, ci1, ci2, fo in zip((self.gen, self.env), (x1, x2),
                                              (y1, y2), (0, 1), (2, 3), (fnt, 12)):
            for n, (x, y) in enumerate(zip(xl, yl)):
                if ((not model[1] and meas['lab'].get() == 1) or
                        (model[1] and not (ci1 == 0 and n not in self.winnerindex))):
                    a.annotate(
                        meas['Name'][n], xy=(x-x3[2], y+y3[3]),
                        fontsize=fo, color=self.color[ci1]
                        )
                if meas['mark'].get() == 1 or model[1]:
                    a.plot(x, y, self.color[ci2], ms=3)
        self.canvas.show()

    def lineplot(self):
        """Plot a linear model for the nominal yield vs environment AMMI2."""
        a = self.a
        a.clear()
        a.set_aspect('auto')
        x_range = self.x_winner[1] - self.x_winner[0]
        y_range = self.y_max - self.y_min
        envmin = min(self.env['IPC1'])
        envmax = max(self.env['IPC1'])
        a.set_xlim(envmin, envmax)
        a.set_ylim(self.y_min-y_range/40, self.y_max+y_range/8)
        a.set_xlabel("Environment IPC1", fontsize=20)
        a.set_ylabel("Nominal Yield", fontsize=20)
        
        for n, (m, i) in enumerate(zip(self.gen['Mean'], self.gen['IPC1'])):
            opt = (self.color[0], 2) if n in self.winnerindex else ('k', .5)
            a.plot([envmin, envmax], [m + i * envmin, m + i *envmax],
                   color=opt[0], linewidth=opt[1])
            
        for n, (lai, x, y) in enumerate(zip(self.winnerindex, self.x_winner, self.y_winner)):
            if n == 0:
                xy = (x + x_range/30, y + y_range/80)
            elif n == 1:
                xy = (x - x_range/7, y + y_range/120)
            else:
                xy = (x + x_range/200, y + y_range/25)
            a.annotate(
                self.gen['Name'][lai], xy=xy, fontsize=17, color=self.color[0]
                )
        self.canvas.show()

if __name__ == '__main__':
    chdir(path.join(getenv('LOCALAPPDATA'), "AMMISOFT"))
    ROOT = tk.Tk()
    Terms(ROOT)
