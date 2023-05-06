import numpy as np
from typing import List,Iterable
from collections import defaultdict
import hashlib 

rev_bitcount_ref = np.array([8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, 7, 6, 6, 5, 6, 5,
       5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4,
       5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 7, 6,
       6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3,
       5, 4, 4, 3, 4, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3,
       3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 7, 6, 6, 5,
       6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4,
       4, 3, 4, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
       5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 6, 5, 5, 4, 5, 4,
       4, 3, 5, 4, 4, 3, 4, 3, 3, 2, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2,
       3, 2, 2, 1, 5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1, 4, 3,
       3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0], dtype=np.int8)

bitmasks = list(range(8))
bitmasks[0] = 0b00000001
bitmasks[1] = 0b00000010
bitmasks[2] = 0b00000100
bitmasks[3] = 0b00001000
bitmasks[4] = 0b00010000
bitmasks[5] = 0b00100000
bitmasks[6] = 0b01000000
bitmasks[7] = 0b10000000

def getbit(bstr,idx):
    byte_idx = idx // 8
    bit_idx = idx % 8
    val = bstr[byte_idx] & bitmasks[bit_idx]
    if val > 0:
        val = 1
    return(val)

def setbit(bstr,idx):
    byte_idx = idx // 8
    bit_idx = idx % 8
    bstr[byte_idx] |= bitmasks[bit_idx]

class DuckTypeError(TypeError):
    pass
                               
class b6hash:
    def __init__(self,input:Iterable[bytes]=None,hashclass=hashlib.blake2b, safety_check:bool=True):
        '''
        Create an instance of the b6hash class to perform similarity 
        hashing of an mset of byte-strings.

        Args:
            input:          An iterable containing byte-strings, representing
                             an mset of byte-strings (default:None) 
            hashclass:      A hash type that conforms with the conventions 
                             of hashlib. (default:hashlib.blake2b)
            safety_check:   If set to True, performs safety checks on input 
                             types which may lower performance. (default:True)
        
        Raises:
            TypeError, ValueError, AttributeError, DuckTypeError
        '''
        self._safety_check = safety_check
        if self._safety_check:
            if not isinstance(input,(Iterable,type(None))):
                raise TypeError(f"'input' must be of type Iterable, not {type(input)}")
            try:
                hashclass().digest_size
                hashclass().update(b'')
                hashclass().digest()
            except AttributeError:
                raise DuckTypeError("'hashclass' must be a type class conforms with the conventions of hashlib. This means instances of hashclass must at least have the 'digest_size' property and functions 'update' and 'digest'")
        self._hashclass = hashclass
        self._OUTSIZE_BYTES = self._hashclass().digest_size
        self._OUTSIZE_VECTOR = self._OUTSIZE_BYTES*8
        
    
        self._vector = np.zeros(self._OUTSIZE_VECTOR,dtype=np.int64)
        self._freq = defaultdict(lambda :0) 
        if not isinstance(input,type(None)):
            self.update(input)

    def update(self,input:Iterable[bytes]) -> np.ndarray:
        '''
        Update this hash object's state with the provided input.

        Args:
            input:  An iterable containing byte-strings, representing
                     an mset of byte-strings 

        Raises:
            TypeError, ValueError
        '''
        if self._safety_check and not isinstance(input,Iterable):
            raise TypeError(f"'input' must be of type Iterable, not {type(input)}. 'input' must be Iterable[bytes]")
        local_freq =  defaultdict(lambda :0)
        for token_b in input:
            if self._safety_check and not isinstance(token_b,bytes):
                raise ValueError("All tokens of 'input' must be of type 'bytes', not {type(token_b)}. 'input' must be Iterable[bytes]")
            local_freq[token_b] += 1     
        for token_b in local_freq:
            previous_count = self._freq[token_b]
            local_count = local_freq[token_b]
            hasher = self._hashclass()
            hasher.update(token_b)
            hasher.update(b'|')
            if previous_count > 0:
                hasher.update(b'1'*previous_count)
            for _ in range(local_count):
                hasher.update(b'1')
                digest = hasher.digest()
                for idx in range(self._OUTSIZE_VECTOR):
                    val = getbit(digest,idx)
                    if val == 1:
                        self._vector[idx] += 1
                    else:
                        self._vector[idx] -= 1
            self._freq[token_b] += local_count

    def vector_digest(self) -> np.ndarray:
        '''
        Get the digest value as an numpy.ndarray object.

        Args:
            n/a
        
        Returns:
            An numpy.ndarray vector of dtype numpy.int64
        '''
        return(self._vector.copy())

    def digest(self) -> bytes:
        '''
        Get the digest value as a bytes object.

        Args:
            n/a
        
        Returns:
            A bytes object representing the digest.
        '''
        output_b = np.zeros(self._OUTSIZE_BYTES,dtype=np.uint8)
        last = 0
        for idx in range(self._OUTSIZE_VECTOR):
            bcnt = self._vector[idx]
            if bcnt > 0:
                setbit(output_b,idx)
            elif bcnt == 0:
                if last == 1:
                    setbit(output_b,idx)
                    last = 0
                else:
                    last = 1
        output_b = output_b.tobytes()
        return(output_b)

    def hexdigest(self) -> str:
        '''
        Get the digest value as a string of hexadecimal digits.

        Args:
            n/a

        Returns:
            A string representing the digest.
        '''
        digest = self.digest()
        s = ''
        for b in digest:
            s += '{:02x}'.format(b)
        return(s)

    @property
    def digest_size(self) -> int:
        '''The number of bytes in a digest produced by .digest()'''
        return(self._OUTSIZE_BYTES)
    
    @property
    def vector_digest_size(self) -> int:
        '''The number of array entries in a vector-digest produced by .hexdigest()'''
        return(self._OUTSIZE_VECTOR)

    @staticmethod
    def bitwise_similarity(b1:bytes,b2:bytes) -> float:
        '''
        Get the bitwise similarity between two byte-strings

        Args:
            b1: A bytes object
            b2: A bytes object
        
        Returns:
            An integer representing the percentage of matching bits.
        '''
        if not isinstance(b1,bytes) or not isinstance(b2,bytes):
            raise TypeError("'b1' and 'b2' must both be of type 'bytes'")
        if len(b1) != len(b2):
            raise ValueError("'b1' and 'b2' must be of the same length")
        count = 0
        for b,bb in zip(b1,b2):
            count += rev_bitcount_ref[b^bb]
        sim = count / (len(b1)*8)
        return(sim)
    
    @staticmethod
    def vector_distance(v1:np.ndarray,v2:np.ndarray,metric:str='euclidean') -> float:
        '''
        Get the distance between two vectors

        Args:
            v1:     An numpy.ndarray object
            v2:     An numpy.ndarray object
            metric: A string representing the distance metric
                     to be used. Either 'euclidean' or 
                     'manhattan' only.
        '''
        if not isinstance(v1,np.ndarray) or not isinstance(v2,np.ndarray):
            raise TypeError("'v1' and 'v2' must both be of type 'numpy.ndarray'")
        if len(v1) != len(v2):
            raise ValueError("'v1' and 'v2' must be of the same length")
        if metric == 'euclidean':
            dist = float(np.sqrt(np.sum(np.power(v1-v2,2))))
        elif metric == 'manhattan':
            dist = float(np.sum(np.abs(v1-v2)))
        else:
            raise ValueError("Currently the only supported methods are 'euclidean' and 'manhattan'")
        return(dist)
    