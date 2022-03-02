import math
from numpy import linalg
from numpy import real


def d_direction(angle, eig1, v1, eig2, v2):
    r1 = abs(math.cos(angle) * v1[0] + math.sin(angle) * v1[1])
    r2 = abs(math.cos(angle) * v2[0] + math.sin(angle) * v2[1])
    return math.sqrt(r1 * r1 * eig1 * eig1 + r2 * r2 * eig2 * eig2)


def add_d(d1, d2, n=4):
    return math.pow((math.pow(d1, -n) + math.pow(d2, -n)) / 2, -1 / n)


def xy_stretch_from_strain(strain):
    return [[strain[0] + 1, strain[1]], [strain[3], strain[4] + 1]]


def stretch_ratio_xy(strain1, strain2, angle_least_compression, angle_most_compression, n=4):
    stretch1 = xy_stretch_from_strain(strain1)
    stretch2 = xy_stretch_from_strain(strain2)
    eigvals, eigvects = linalg.eig(stretch1)
    eig1 = eigvals[0]
    eig2 = eigvals[1]
    m = max(eig1, eig2)  # normalize longer eigenvector to 1
    eig1 = eig1 / m
    eig2 = eig2 / m
    v1 = eigvects[0]
    v2 = eigvects[1]
    eigvals, eigvects = linalg.eig(stretch2)
    eig3 = eigvals[0]
    eig4 = eigvals[1]
    m = max(eig3, eig4)  # normalize longer eigenvector to 1
    eig3 = eig3 / m
    eig4 = eig4 / m
    v3 = eigvects[0]
    v4 = eigvects[1]
    dmax = add_d(d_direction(angle_least_compression, eig1, v1, eig2, v2),
                 d_direction(angle_least_compression, eig3, v3, eig4, v4), n=n)
    dmin = add_d(d_direction(angle_most_compression, eig1, v1, eig2, v2),
                 d_direction(angle_most_compression, eig3, v3, eig4, v4), n=n)
    return dmax / dmin


def least_compressive_direction_xy(strain1, strain2, n=4):
    stretch_xy_1 = xy_stretch_from_strain(strain1)
    stretch_xy_2 = xy_stretch_from_strain(strain2)
    eigvals, eigvects = linalg.eig(stretch_xy_1)
    eig1 = eigvals[0]
    eig2 = eigvals[1]
    m = max(eig1, eig2)  # normalize longer eigenvector to 1
    eig1 = eig1 / m
    eig2 = eig2 / m
    v1 = eigvects[0]
    v2 = eigvects[1]
    eigvals, eigvects = linalg.eig(stretch_xy_2)
    eig3 = eigvals[0]
    eig4 = eigvals[1]
    m = max(eig3, eig4)  # normalize longer eigenvector to 1
    eig3 = eig3 / m
    eig4 = eig4 / m
    v3 = eigvects[0]
    v4 = eigvects[1]
    # First, find grossly the minimum by a first gross screen
    angles_to_screen = []
    for a_index in range(901):
        angles_to_screen.append(a_index / 900 * math.pi / 2)
    angle = angles_to_screen[0]
    max_d = add_d(d_direction(angle, eig1, v1, eig2, v2), d_direction(angle, eig3, v3, eig4, v4), n=n)
    max_angle = angle
    for theAngle in angles_to_screen:
        new_d = add_d(d_direction(theAngle, eig1, v1, eig2, v2), d_direction(theAngle, eig3, v3, eig4, v4), n=n)
        if new_d > max_d:
            max_angle = theAngle
            max_d = new_d
    return max_angle


def most_compressive_direction_xy(strain1, strain2, n=4):
    stretch_xy_1 = xy_stretch_from_strain(strain1)
    stretch_xy_2 = xy_stretch_from_strain(strain2)
    eigvals, eigvects = linalg.eig(stretch_xy_1)
    eig1 = eigvals[0]
    eig2 = eigvals[1]
    m = max(eig1, eig2)  # normalize longer eigenvector to 1
    eig1 = eig1 / m
    eig2 = eig2 / m
    v1 = eigvects[0]
    v2 = eigvects[1]
    eigvals, eigvects = linalg.eig(stretch_xy_2)
    eig3 = eigvals[0]
    eig4 = eigvals[1]
    m = max(eig3, eig4)  # normalize longer eigenvector to 1
    eig3 = eig3 / m
    eig4 = eig4 / m
    v3 = eigvects[0]
    v4 = eigvects[1]
    # First, find grossly the minimum by a first gross screen
    angles_to_screen = []
    for a_index in range(901):
        angles_to_screen.append(a_index / 900 * math.pi / 2)
    angle = angles_to_screen[0]
    min_d = add_d(d_direction(angle, eig1, v1, eig2, v2), d_direction(angle, eig3, v3, eig4, v4), n=n)
    min_angle = angle
    for theAngle in angles_to_screen:
        new_d = add_d(d_direction(theAngle, eig1, v1, eig2, v2), d_direction(theAngle, eig3, v3, eig4, v4), n=n)
        if new_d < min_d:
            min_angle = theAngle
            min_d = new_d
    return min_angle
