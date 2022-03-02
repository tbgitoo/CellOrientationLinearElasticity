import math
from numpy import linalg
from numpy import real




# These functions need strain values, the identity contribution is added here
def eigenvalue_ratio(s):
    e=linalg.eigvals([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]], [s[6], s[7], 1 + s[8]]])
    return max(real(e))/min(real(e))

def eigenvalue_ratio_xy(s):
    e=linalg.eigvals([[1 + s[0], s[1]], [s[3], 1 + s[4]]])
    return max(real(e))/min(real(e))



def most_compressive_eigenvector(s):
    e=linalg.eigvals([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])
    if(e[0]<e[1] and e[0]<e[2]):
        return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])[1][0]
    else:
        if(e[0]<e[1]):
            return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])[1][1]
        else:
            return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]], [s[6], s[7], s[8] + 1]])[1][1]

def least_compressive_eigenvector(s):
    e=linalg.eigvals([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])
    if(e[0]>e[1] and e[0]>e[2]):
        return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])[1][0]
    else:
        if(e[0]>e[1]):
            return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]],[s[6],s[7],s[8]+1]])[1][1]
        else:
            return linalg.eig([[1 + s[0], s[1], s[2]], [s[3], 1 + s[4], s[5]], [s[6], s[7], s[8] + 1]])[1][1]


def most_compressive_eigenvector_xy(s):
    e=linalg.eigvals([[1 + s[0], s[1]], [s[3], 1 + s[4]]])
    if(e[0]<e[1]):
        return linalg.eig([[1 + s[0], s[1]], [s[3], 1 + s[4]]])[1][0]
    else:
        return linalg.eig([[1 + s[0], s[1]], [s[3], 1 + s[4]]])[1][1]

def least_compressive_eigenvector_xy(s):
    e=linalg.eigvals([[1 + s[0], s[1]], [s[3], 1 + s[4]]])
    if(e[0]>e[1]):
        return linalg.eig([[1 + s[0], s[1]], [s[3], 1 + s[4]]])[1][0]
    else:
        return linalg.eig([[1 + s[0], s[1]], [s[3], 1 + s[4]]])[1][1]









