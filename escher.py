from cmath import sin, cos
#import scipy.misc as spmisc
import PIL.Image

epsilon = 1e-10
pi = 3.14159265359

def iround(x):
    return int(x + 0.5)

class MoebiusTransformation(object):
    def __init__(self, z_0, z_1):
        self.z_0 = z_0
        z_1_image= (z_1 - z_0) / (1.0 - z_0.conjugate() * z_1)
        self.eps = abs(z_1_image) / z_1_image

    def __call__(self, z):
        return self.eps * (z - self.z_0) / (1.0 - self.z_0.conjugate() * z)

    def inverse(self):
        z_0 = self(0.0)
        z_1 = self(0.5)
        return MoebiusTransformation(z_0, z_1)

class MarkedHalfPlane(object):
    def __init__(self, z_1, z_2, z):
        if abs(z_1 - z_2) < epsilon:
            raise ValueError("A line should be defined by two different "
                    "points.")
        self.z_1 = z_1
        self.z_2 = z_2
        self.z = z
        self.moebius = MoebiusTransformation(z_1, z_2)
        self.moebius_inv = self.moebius.inverse()
        self.w_1 = self.moebius(z_1)
        self.w_2 = self.moebius(z_2)
        self.w = self.moebius(z)
        if self.w.imag < epsilon:
            raise ValueError("The marking point of the MarkedHalfplane should not "
                    "be on the boundary line.")

    def reflect_in(self, z):
        s = self.moebius(z)
        if s.imag * self.w.imag >= 0:
            return z
        return self.moebius_inv(s.conjugate())

    def sines_of_angles(self, z):
        s = self.moebius(z)
        ratios = (self.w / s, s)
        return (x.imag/ abs(x) for x in ratios)

def get_coordinates(halfplanes, z):

    #if z is a vertex, sines_of_angles would lead to division by zero
    for i in range(3):
        if abs(halfplanes[i].z_1 - z) < epsilon:
            l = [0.0, 0.0, 0.0]
            l[i] = 1.0
            return l

    weights = [[1.0]*3 for i in range(3)]
    for i in range(3):
        (a, b) = halfplanes[i].sines_of_angles(z)
        weights[(i+1)%3][ (i+2)%3] = b / (a + b)
        weights[(i+2)%3][ (i+1)%3] = a / (a + b)

    #find vertex with largest weight (to avoid dividing by zero)
    for i in range(3):
        if all(weights[i][(i+j)%3] < weights[(i+j)%3][i] + epsilon for 
                j in {1, 2}):
            lw = i
        
    coord = [weights[lw][i]/weights[i][lw] for i in range(3)]
    return [x/sum(coord) for x in coord]


def linear_combination(weights, vertices):
    return (iround(sum(weights[i] * vertices[i][j] for i in range(3))) 
            for j in range(2))

def complex_to_pixels(z, halfsize):
    return (iround(halfsize * (1.0 - z.imag)),iround(halfsize * (1.0 + z.real)))

def pixels_to_complex(pair, halfsize):
    return complex(pair[1]/halfsize - 1.0, 1.0 - pair[0]/halfsize)
    

white = (255, 255, 255)

class ImageSource(object):
    def __init__(self, image_file = None, vertices = None, default_color = white):
        self.default_color = default_color
        if image_file is None:
            self.array = PIL.Image.new('RGB', (1, 1), default_color).load()
            self.vertices = ((0,0),(0,0),(0,0))
        else:
            self.vertices = vertices
            #self.array = spmisc.imread(image_file).tolist()
            self.array = PIL.Image.open(image_file).convert('RGB').load()


def create_image(output_file_name, size, p, q, r, image_sources):
    n = len(image_sources)
    if not all((x % n) == 0 for x in {p, q, r}):
        raise ValueError("The number of triangles around each vertex "
                "must be divisible by the number of image sources. ")
    vertices = triangle(p, q, r)
    new_image = [[[0]*3 for i in range(size)] for j in range(size)]
    temp_images = [[[[-1]*3 for i in range(size)] for j in range(size)] for x in image_sources]
    halfplanes = [MarkedHalfPlane(vertices[i], vertices[(i+1)%3], vertices[(i+2)%3])
            for i in range(3)]
    
    halfsize = size/2.0
    for i in range(size):
        for j in range(size):
            z = pixels_to_complex((i, j), halfsize)
            
            if abs(z) >= 1:
                new_image[i][j] = white
                continue

            do_more = True
            count = 0
            while do_more:
                do_more = False
                for k in range(3):
                    newz = halfplanes[k].reflect_in(z)
                    if newz != z:
                        do_more = True
                        z = newz
                        count += 1
            
            (newi, newj) = complex_to_pixels(z, halfsize)
            current_temp_image = temp_images[count % n]
            current_image_source = image_sources[count % n]
            #print current_temp_image[newi][newj]
            if current_temp_image[newi][newj][0] < 0: #pixel is uninitialized
                weights = get_coordinates(halfplanes, z)
                (x, y) = linear_combination(weights, image_sources[count % n].vertices)
                try:
                    current_temp_image[newi][newj] = current_image_source.array[x,y]
                except IndexError: #if the coordinates are out of range
                    current_temp_image[newi][newj] = current_image_source.default_color
            new_image[i][j] = current_temp_image[newi][newj]

    #spmisc.imsave(output_file_name, new_image)
    new_img = PIL.Image.new("RGB", (size, size))
    new_img.putdata([x for row in new_image for x in row])
    new_img.save(output_file_name)

def cosh_side_length(alpha, beta, gamma):
    return (cos(alpha)*cos(beta) + cos(gamma)) / (sin(alpha) * sin(beta))

def point_from_origin(cosh_of_dist):
    return ((cosh_of_dist - 1.0) / (cosh_of_dist + 1.0))**(0.5)

def triangle(p, q, r):
    if 1.0/p + 1.0/q + 1.0/r >= 0.5 - epsilon:
        raise ValueError("The sum of angles in a hyperbolic triangle is less"
                " than pi, so 1/p + 1/q + 1/r < 1/2 should be satisfied.")
    alpha = 2*pi/p
    beta = 2*pi/q
    gamma = 2*pi/r

    return (0.0, point_from_origin(cosh_side_length(alpha, beta, gamma)),
            point_from_origin(cosh_side_length(alpha, gamma, beta)) *
            complex(cos(alpha), sin(alpha)))

image_source = ImageSource("stinkbug.png", ((0,800), (0, 0), (800, 0)))
image_source2 = ImageSource(default_color = (50, 223, 146))
create_image("circle.png", 1000, 6, 6, 8, [image_source, image_source2])



