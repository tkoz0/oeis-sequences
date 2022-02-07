'''
TPTreePreOrder
- generator object for (length,number) in pre order
TPTreePostOrder
- generator object for (length,number,num_children) in post order
TPTreeFull
- returns an object representing the whole tree (numbers only)

Constructor arguments:
__init__(self, ptype, base, root, buf)

Call .generator() to get the object it is setup for
'''

from typing import Any, BinaryIO, Callable, Generator, List, Tuple
import sys

# (length, prime) (for pre order, number of child nodes not known yet)
TreeNodePre = Tuple[int,int]

# (length, prime, children) (post order, after number of children is known)
TreeNodePost = Tuple[int,int,int]

# numbers only
TPTree = Tuple[int,List['TPTree']]

# maps (length,number) with tree value bytes to another (length,number)
# length is digits in the base
NextPrime = Callable[[int,int,bytes],Tuple[int,int]]

# end byte
END = b'\xff'

# reads a single byte, erroring if unable
def _read_next_byte(buf: BinaryIO = sys.stdin.buffer) -> bytes:
    b = buf.read(1)
    assert b
    return b

# returns b'\xff' for end byte, otherwise size bytes for the value
def _read_next_value(size: int, buf: BinaryIO = sys.stdin.buffer) -> bytes:
    b = _read_next_byte(buf)
    if b == END:
        return b
    for i in range(1,size):
        b += _read_next_byte(buf)
    return b

# counts digits of v in base b
def _count_digits(v: int, b: int) -> int:
    assert v >= 0 and b >= 2
    d = 0
    while v:
        v //= b
        d += 1
    return d

class TPTreeBase:
    ptype: str
    base: int
    root: int
    buf: BinaryIO
    vsize: int
    nextp: NextPrime
    def __init__(self, ptype: str, base: int, root: int = 0,
            buf: BinaryIO = sys.stdin.buffer):
        assert ptype in ['r','l','lor','lar']
        self.ptype = ptype
        assert 2 <= base < 256
        self.base = base
        assert root >= 0
        self.root = root
        self.buf = buf
        self.vsize = 1+(len(ptype)//2) # 1 for r,l and 2 for lor,lar
        self.nextp = \
        {
            'r': self._next_r,
            'l': self._next_l,
            'lor': self._next_lor,
            'lar': self._next_lar
        }[ptype]
    def _next_r(self, l: int, v: int, b: bytes) -> Tuple[int,int]:
        assert 0 < b[0] < self.base
        return (l+1,self.base*v+b[0])
    def _next_l(self, l: int, v: int, b: bytes) -> Tuple[int,int]:
        assert 0 < b[0] < self.base
        return (l+1,v+(10**l)*b[0])
    def _next_lor(self, l: int, v: int, b: bytes) -> Tuple[int,int]:
        assert 0 < b[1] < self.base
        if b[0] == 0:
            return self._next_l(l,v,b[1:])
        elif b[0] == 1:
            return self._next_r(l,v,b[1:])
        else:
            assert 0, 'invalid lor first byte'
            return 0,0
    def _next_lar(self, l: int, v: int, b: bytes) -> Tuple[int,int]:
        if l == 0: # root, allows single digit primes
            assert b != '\x00\x00', 'lar digits both zero'
            return (2 if b[0] else 1,self.base*b[0]+b[1])
        else:
            assert 0 < b[0] < self.base and 0 < b[1] < self.base
            return (l+2,self.base*v+b[1]+(10**(l+1))*b[0])
    def generator(self) -> Any:
        assert 0, 'must override generator'

# generator object for (length,number) in pre order
class TPTreePreOrder(TPTreeBase):
    def _generator(self, l: int, n: int) -> Generator[TreeNodePre,None,None]:
        yield (l,n)
        while True:
            v = _read_next_value(self.vsize,self.buf)
            if v == END:
                break
            l2,n2 = self.nextp(l,n,v)
            yield from self._generator(l2,n2)
    def generator(self) -> Generator[TreeNodePre,None,None]:
        for _ in range(self.vsize): # skip root
            _read_next_byte(self.buf)
        yield from self._generator(_count_digits(self.root,self.base),self.root)

# generator object for (length,number,children) in post order
class TPTreePostOrder(TPTreeBase):
    def _generator(self, l: int, n: int) -> Generator[TreeNodePost,None,None]:
        c = 0
        while True:
            v = _read_next_value(self.vsize,self.buf)
            if v == END:
                break
            c += 1
            l2,n2 = self.nextp(l,n,v)
            yield from self._generator(l2,n2)
        yield (l,n,c)
    def generator(self) -> Generator[TreeNodePost,None,None]:
        for _ in range(self.vsize): # skip root
            _read_next_byte(self.buf)
        yield from self._generator(_count_digits(self.root,self.base),self.root)

# creates the full tree object with values only
class TPTreeFull(TPTreeBase):
    def _generator(self, l: int, n: int) -> TPTree:
        children: list[TPTree] = []
        while True:
            v = _read_next_value(self.vsize,self.buf)
            if v == END:
                break
            l2,n2 = self.nextp(l,n,v)
            children.append(self._generator(l2,n2))
        return (n,children)
    def generator(self) -> TPTree:
        for _ in range(self.vsize): # skip root
            _read_next_byte(self.buf)
        return self._generator(_count_digits(self.root,self.base),self.root)
