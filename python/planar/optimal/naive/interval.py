from math import *

class Interval:
    """Represents real intervals."""

    def __init__(self, left, right, left_closed=False, right_closed=False):
        if left < right:
            self.left = left
            self.right = right
            self.left_closed = left_closed
            self.right_closed = right_closed
        elif left > right:
            self.left = None
            self.right = None
            self.left_closed = None
            self.right_closed = None
        elif left_closed and right_closed:
            self.left = left
            self.right = right
            self.left_closed = left_closed
            self.right_closed = right_closed
        else:
            self.left = None
            self.right = None
            self.left_closed = None
            self.right_closed = None

    """The inequality operators establish a partial ordering based on
       strict inequality. == and != don't follow the same rules. a >= b
       and a <= b implies that the intervals overlap, but does not imply
       a == b."""

    def __lt__(self, other):
        if self.is_empty() or other.is_empty():
            return False
        else:
            return (self.right < other.left or
                   (self.right == other.left and
                   not (self.right_closed or other.left_closed)))

    def __le__(self, other):
        return not (self > other)

    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return not (self < other)

    def __key(self):
        return (self.left, self.right, self.left_closed, self.right_closed)

    def __eq__(self, other):
        return type(self) == type(other) and self.__key() == other.__key()

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.__key())

    def __str__(self):
        if self.is_empty():
            return "(0, 0)"
        build = ""
        if self.left_closed:
            build += '['
        else:
            build += '('
        build += str(self.left)
        build += ", "
        build += str(self.right)
        if self.right_closed:
            build += ']'
        else:
            build += ')'
        return build

    def __repr__(self):
        return str(self)

    def is_empty(self):
        return self == EMPTY_INTERVAL

    def intersection(self, other):
        """Returns an intersection with another interval, or None if the
           the intersection is empty."""
        if self.is_empty() or other.is_empty():
            return EMPTY_INTERVAL
        if self < other or other < self:
            return EMPTY_INTERVAL
        if self.left < other.left:
            left = other.left
            left_closed = other.left_closed
        elif self.left > other.left:
            left = self.left
            left_closed = self.left_closed
        else:
            left = self.left
            left_closed = self.left_closed and other.left_closed
        if self.right > other.right:
            right = other.right
            right_closed = other.right_closed
        elif self.right < other.right:
            right = self.right
            right_closed = self.right_closed
        else:
            right = self.right
            right_closed = self.right_closed and other.right_closed
        return Interval(left, right, left_closed, right_closed)

    def union(self, other):
        """Returns a union with another interval, or None if the
           intervals are disjoint (use IntervalUnion for unions of
           disjoint intervals)"""
        if self.is_empty():
            return other
        if other.is_empty():
            return self
        if self < other or other < self:
            return None
        if self.left < other.left:
            left = self.left
            left_closed = self.left_closed
        elif self.left > other.left:
            left = other.left
            left_closed = other.left_closed
        else:
            left = self.left
            left_closed = self.left_closed or other.left_closed
        if self.right > other.right:
            right = self.right
            right_closed = self.right_closed
        elif self.right < other.right:
            right = other.right
            right_closed = other.right_closed
        else:
            right = self.right
            right_closed = self.right_closed or other.right_closed
        return Interval(left, right, left_closed, right_closed)

EMPTY_INTERVAL = Interval(0, 0)

class IntervalUnion:
    """Represents a union of real intervals. Not well optimized,
       refactor code if many disjoint unions are needed."""
    def __init__(self, intervals=[]):
        self.intervals = []
        for interval in intervals:
            self.unionadd(interval)

    def __eq__(self, other):
        return self.intervals == other.intervals

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return str(self.intervals)

    def __repr__(self):
        return str(self)

    def unionadd(self, new_interval):
        if new_interval.is_empty():
            return
        result = []
        for interval in self.intervals:
            if new_interval is None:
                result.append(interval)
            else:
                if interval < new_interval:
                    result.append(interval)
                elif new_interval < interval:
                    result.append(new_interval)
                    new_interval = None
                    result.append(interval)
                else:
                    new_interval = interval.union(new_interval)
        if new_interval is not None:
            result.append(new_interval)
        self.intervals = result

    def union(self, other):
        new_union = IntervalUnion(self.intervals)
        for interval in other.intervals:
            new_union.unionadd(interval)
        return new_union

TAU = 2*pi

def make_angle_range(start, end, start_closed=False, end_closed=False):
    """Creates an IntervalUnion representing a set of angles.
       Range of values: [0, 2pi)"""
    start = start % TAU
    end = end % TAU
    if start < end:
        return IntervalUnion([Interval(start, end, left_closed=start_closed, right_closed=end_closed)])
    elif end < start:
        return IntervalUnion([Interval(0, end, left_closed=True, right_closed=end_closed), Interval(start, TAU, left_closed=start_closed, right_closed=False)])
    else:
        return IntervalUnion()

FULL_ANGLE_RANGE = make_angle_range(0, TAU, start_closed=True)

#TEST

a = Interval(0, 10)
b = Interval(-5, 5)
c = a.union(b)
