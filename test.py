import matplotlib.pyplot as plt
import matplotlib

class Interactivity:
    def __init__(self, fig = None):
        self.fig = plt.gcf() if fig is None else fig
        self.ax = self.fig.gca()
        self.connections = ()

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()

    def connect(self):
        """ Install the event handlers for the plot. """
        self.connections = (
            self.fig.canvas.mpl_connect('button_press_event', self.onclick),
            self.fig.canvas.mpl_connect('pick_event', self.onpick),
            self.fig.canvas.mpl_connect('key_press_event', self.on_key),
        )

    def disconnect(self):
        """ Uninstall the event handlers for the plot. """
        for connection in self.connections:
            self.fig.canvas.mpl_disconnect(connection)

    def draw_line(self, startx, starty):
        xy = plt.ginput(1)
        x = [startx, xy[0][0]]
        y = [starty, xy[0][1]]
        self.ax.plot(x, y, picker=True, pickradius=5, color='blue')
        self.ax.figure.canvas.draw_idle()

    def onclick(self, event):
        """
        This implements click functionality. If it's a double click do
        something, else ignore.
        Once in the double click block, if its a left click, wait for a further
        click and draw a line between the double click co-ordinates and that
        click (using ginput(1) - the 1 means wait for one mouse input - a
        higher number is used to get multiple clicks to define a polyline)
        """
        print('onclick')
        if event.dblclick:
            if event.button == 1:
                print('event location', event.xdata, event.ydata)
                self.draw_line(event.xdata, event.ydata)


    def onpick(self, event):
        """
        Handles the pick event - if an object has been picked, store a
        reference to it.  We do this by simply adding a reference to it
        named 'picked_object' to the axes object.
        """
        print('onpick')
        this_artist = event.artist
        print('this_artist', this_artist)
        # the picked object is available as event.artist
        self.ax.picked_object = this_artist
        print('picked_object', self.ax.picked_object)

    def on_key(self, event):
        """
        Function to be bound to the key press event
        If the key pressed is delete and there is a picked object,
        remove that object from the canvas
        """
        print('onkey: ', event.key)
        if event.key == 'delete' and self.ax.picked_object:
            self.ax.picked_object.remove()
            self.ax.picked_object = None
            self.ax.figure.canvas.draw_idle()
            

# Basic plot to test the functionality
fig = plt.figure(figsize = (10,10))
gs = fig.add_gridspec(10,10)
ax101 = fig.add_subplot(gs[:,:])
ax101.set_ylim(0,10)
ax101.set_xlim(0,10)
with Interactivity():
    plt.show()