"""
window.py: define the basic class to describe the windows.
"""
import pickle


class Window(object):
    """
    Class defines a single window used
    """

    def __init__(self, left=None, right=None, channel=None, network=None, gcmtid=None, station=None, phases=None):
        super().__init__()
        self.left = left
        self.right = right
        self.channel = channel
        self.network = network
        self.gcmtid = gcmtid
        self.station = station
        self.phases = phases

    def __repr__(self):
        return f"Windows(left={self.left},right={self.right},channel={self.channel},network={self.network},gcmtid={self.gcmtid},station={self.station},phases={self.phases})"

    def __add__(self, others):
        left = min(self.left, others.left)
        right = max(self.right, others.right)
        assert self.channel == others.channel
        assert self.network == others.network
        assert self.station == others.station
        assert self.gcmtid == others.gcmtid
        self.phases = self.phases + others.phases
        return Window(left=left, right=right, channel=self.channel, network=self.network, station=self.station, phases=self.phases+others.phases, gcmtid=self.gcmtid)


class Windows_collection(object):
    def __init__(self, newwindows=None):
        if(newwindows == None):
            self.windows = []
        else:
            self.windows = newwindows.windows

    def append_window(self, window):
        self.windows.append(window)

    def save(self, fname):
        # save windows
        with open(fname, "wb") as handle:
            pickle.dump(self.windows, handle)

    def load(self, fname):
        with open(fname, "rb") as handle:
            self.windows = pickle.load(handle)

    def combine(self, others):
        newwindows = self.windows + others.windows
        return Windows_collection(newwindows=newwindows)

    def merge_windows(self):
        # merge windows
        self.windows = sorted(self.windows, key=lambda item: item.left)
        new_windows = [self.windows[0]]
        for index in range(1, len(self.windows)):
            current_window = new_windows[-1]
            if (current_window.right > self.windows[index].left):
                new_windows[-1] = new_windows[-1]+self.windows[index]
            else:
                new_windows.append(self.windows[index])
        self.windows = new_windows
