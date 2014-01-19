"""
Written by: Balazs Strenner, UW-Madison, Mathematics Department
Initial version in C++ in 2006.
Current version with improved algorithm: 2014 January.

This script is for creating Escher-like tilings of the Poincare disk.

MATHEMATICAL BACKGROUND:

    The Poincare disk is a model of hyperbolic geometry, where the points
    are the points of the unit disk and lines are either lines containing
    the origin, or circlular arcs intersecting the unit circle in right 
    angle. A hyperbolic triangle is therefore a region bounded by three 
    line segments (that either look straight or circular) obtained from
    three points. The angle between two hyperbolic lines is the same as the
    usual Euclidean angle. A fun fact is that the sum of the angles of a
    hyperbolic triangle is always less than 180 degrees, and this sum can
    take any value between 0 and 180 degrees.

    On the Euclidean plane, if one takes an equilateral triangle and keeps
    reflecting it to different sides of the triangle, one can obtain a
    tiling of the whole plane using the same triangle (without overlapping
    and leftover uncovered areas). One can do this with a 90-60-30, and
    a 120-30-30 triangle as well, but that's all. Indeed, each angle of the
    triangle must be a quotient of 360 degrees by an positive integer,
    otherwise the triangles overlap when flipping them around a vertex.
    Since the sum of the three angles is 180 degrees, we have the equation
    360/a + 360/b + 360/c = 180 where a, b, c are positive integers, at least
    3. I.e. 1/a + 1/b + 1/c = 1/2. One can easily check that only the triples
    (3, 6, 6), (4, 6, 12) and (6, 6, 6) work as (a, b, c).

    On the hyperbolic plane however, the required condition is
    1/a + 1/b + 1/c < 1/2, and here there are infinitely many possibilities,
    which explains why the tessallations of the hyperbolic plane are much 
    richer than that of the Euclidean plane. 

    If all of a, b, c are even, then the tilings have the nice property that 
    the any side of any triangle in the tessallation is a subset of a line
    that run along the boundaries of triangles. If any of a, b, c is odd,
    then some of these lines cut into the interior of other triangles which
    makes the situation more complicated. Therefore we restrict ourselves
    to the case when (a, b, c) = (2p, 2q, 2r) for some integers p, q, r.
    In the hyperbolic plane, 1/p + 1/q + 1/r < 1 should hold.

WHAT THE PROGRAM DOES:

    Creates a tiling of the Poincare disk using any triangle shape that 
    satisfies 1/p + 1/q + 1/r < 1. The inside of each hyperbolic triangle
    is filled with an image. This image should be specified by the name
    of an image file and the coordinates of the three vertices of the
    triangular region that is then pasted into the hyperbolic triangle.
    One can actually use multiple images. If one uses two images, they
    will fill the triangles alternatingly. (E.g. one image with Tom,
    another with Jerry, and the infinitely Toms will chase infinitely many
    Jerrys in the Poincare disk.)If more image, however such
    guarantees are impossible to pose, so the outcome will look a little
    random. Not even that is guaranteed that no two neighboring triangle
    will be filled using the same image.

    When specifying the image file and the triangle, the vertices of the 
    triangle are allowed to be outside of the image, and in this case any
    pixel inside the bounded triangle that is outside of the image will
    get a default color (that can also be specified). This is useful if
    there is little margin around, say, a face on an image, and thus it is
    hard to fully enclose it in a triangle that is contained in the image.
    
    One can also omit specifying an image file. In that case the triangles
    will be colored to the default color.

REQUIREMENTS:

    - Python Image Library (PIL): since it is not maintained anymore,
        check out its fork "pillow". (If you have Sage, you already have
        PIL installed.)

EXAMPLE USE:

    Have a stinkbug.png image of size 800x600 in your current directory. 
    Open your favorite
    Python or Sage interpreter (python, ipython, pypy, sage, etc. in a 
    UNIX terminal), then type:

        >>> import escher 
        >>> triangle_image = escher.TriangleImage("stinkbug.png", ((0,599), (0, 0), (799, 0)))
        >>> triangle_image2 = escher.TriangleImage(default_color = (50, 223, 146))
        >>> escher.create_image([triangle_image, triangle_image2], 3, 3, 4, 500, "output.png") # doctest: +SKIP

    It is possible to use a rectangular image for the tessallation. Simply
    create two TriangleImage objects from the same image, one from the 
    upper left, another from the lower right triangle:

        >>> escher.create_image_from_rectangle("stinkbug.png", 5, 5, 5, 500, "output.png") # doctest: +SKIP

    The above command is a shorthand for the following three together which use the
    upper left and bottom right triangles for the image to simulate a rectange:

        >>> triangle_image = escher.TriangleImage("stinkbug.png", ((0,599), (0, 0), (799, 0)))
        >>> triangle_image2 = escher.TriangleImage("stinkbug.png", ((0,599), (799, 599), (799, 0)))
        >>> escher.create_image([triangle_image, triangle_image2], 3, 3, 4, 500, "output.png") # doctest: +SKIP

    To get more help about how to define TriangleImage objects and how to
    use the create_image function, type:

        >>> TriangleImage? # doctest: +SKIP

    or

        >>> create_image? # doctest: +SKIP

"""







from math import sin, cos, pi
import PIL.Image

_epsilon = 1e-10

def _iround(x):
    """ Rounds a float to the nearing integer.  """
    return int(x + 0.5)

class MoebiusTransformation(object):
    """
    A Moebius transformation of fixing the unit disk in the complex plane.

    All such transformations have the form 
    f(z) = eps*(z-z_0)/(1-z_0.conjugate()*z) with some |z_0| < 1 and 
    |eps| = 1. Moebius transformations are exactly the orientation-
    preserving isomeries of the hyperbolic plane as the Poincare disk.

    INPUT:

    ``z_0`` - The z_0 as above, i.e. the complex number mapped into 0.

    ``z_1`` - A complex number in the unit disk, different from z_0, that
        will be mapped to a positive real number. This parameter uniquely
        determines ``eps`` above.

    EXAMPLES:

    The following is the identity map, so its inverse is the identity, too:
        
        >>> f = MoebiusTransformation(0, 0.5); f(-0.3)
        -0.3
        >>> g = f.inverse(); g(-0.3)
        -0.3

    The following is a rotation about the origin by 90 degrees in the 
    counterclockwise direction. Its inverse is the rotation in the other
    direction:

        >>> f = MoebiusTransformation(0, complex(0,-0.5)); f(-0.3)
        -0.3j
        >>> g = f.inverse(); g(-0.3)
        0.3j

    The following is a more complex one:

        >>> f = MoebiusTransformation(-0.2, complex(0.1, 0.2)); f(-0.3)
        (-0.09076037023968804+0.05549678689815319j)
        >>> g = f.inverse(); g(-0.3)
        (-0.43778641728160506-0.13584442237179067j)

    """
    def __init__(self, z_0, z_1):
        self.z_0 = z_0
        z_1_image= (z_1 - z_0) / (1.0 - z_0.conjugate() * z_1)
        self.eps = abs(z_1_image) / z_1_image

    def __call__(self, z):
        """ Calculates the image of ``z`` under the transformation self.  """
        return self.eps * (z - self.z_0) / (1.0 - self.z_0.conjugate() * z)

    def inverse(self):
        """ Returns the inverse transformation of self.  """
        z_0 = self(0.0)
        z_1 = self(0.5)
        return MoebiusTransformation(z_0, z_1)



class _MarkedHalfPlane(object):
    """
    A halfplane in the hyperbolic plane marked with two points on
    its boundary and one point inside.

    Defining a line, a halfplane and a triangle class separately would
    surely be more useful for reusability in a library, but when
    generating a large image, defining the same Moebius transformation
    over and over is costly. When generating the image, we reflect different
    points through the same three lines over and over (there is also the
    calculation of sines of angles), thus by this
    more complex class we get a faster code.

    INPUT:

    ``z_1``, ``z_2`` - the two complex numbers on the boundary of the 
        halfplane

    ``z_3`` - the complex number inside the halfplane

    EXAMPLES:

    This is the halfplane corresponding to complex numbers with 
    non-negative imaginary part:

        >>> m = _MarkedHalfPlane(0, 0.5, complex(0, 0.5))

    """
    def __init__(self, z_1, z_2, z_3):
        if abs(z_1 - z_2) < _epsilon:
            raise ValueError("A line should be defined by two different "
                    "points.")
        self.z_1 = z_1
        self.z_2 = z_2
        self.moebius = MoebiusTransformation(z_1, z_2)
        self.moebius_inv = self.moebius.inverse()
        self.w = self.moebius(z_3)
        if self.w.imag < _epsilon:
            raise ValueError("The marking point of the _MarkedHalfplane should not "
                    "be on the boundary line.")

    def reflect_in(self, z):
        """
        Reflects a point into the halfplane through its boundary line.

        INPUT:

        ``z`` - a complex number in the Poincare disk

        OUTPUT:

        - a complex number -

        EXAMPLES:

        If the point is already in the halfplane, the output is the same
        point:

            >>> f = _MarkedHalfPlane(0, 0.5, complex(0, 0.5))
            >>> f.reflect_in(-0.3)
            -0.3
            >>> f.reflect_in(complex(-0.2, 0.1))
            (-0.2+0.1j)

        In the currect case, since the boundary line lies on the real axis,
        reflecting is the same as conjugation:

            >>> f.reflect_in(complex(0.3, -0.1))
            (0.3+0.1j)

        """
        s = self.moebius(z)
        if s.imag * self.w.imag > -_epsilon:
            return z
        return self.moebius_inv(s.conjugate())

    def sines_of_angles(self, z):
        """
        Returns the tuple of the sines of the angles (z_3 z_1 z) and
        (z_2 z_1 z).

        This function is used by _get_coordinates to calculate the 
        weight-coordinates of a point with respect to z_1, z_2, z_3.

        INPUT:

        - ``z`` - a complex number

        OUTPUT:

        - a tuple of two real numbers

        EXAMPLES:

        The sine of 45 degrees is sqrt(2)/2 = 0.7071067811865475, so:

            >>> m = _MarkedHalfPlane(0, 0.5, complex(0, 0.5))
            >>> [round(x, 10) for x in m.sines_of_angles(complex(0.1, 0.1))]
            [0.7071067812, 0.7071067812]
        """
        s = self.moebius(z)
        ratios = (self.w / s, s)
        return tuple(x.imag/ abs(x) for x in ratios)

def _get_coordinates(halfplanes, z):
    """
    Returns the weight-coordinates of a point with respect to the
    vertices of a triangle.

    Here is the motivation. Given three points x, y, z on the Euclidean
    plane, one can parametrize the points inside a triangle by triples
    (a, b, c) of non-negative real numbers with a + b + c = 1, since
    the linear combination a*x + b*y + c*z defines a point inside the
    triangle. If the triangle is non-degenerate, this parametrization is
    unique.

    Linear combinations don't work on the hyperbolic plane unfortunately,
    but Ceva's theorem does. So there is a straightfoward way to get the 
    parameters (a, b, c) by looking at the ratios the three angle bisectors
    divide the opposite sides in, and this works in in hyperbolic 
    geometry, too.

    INPUT:

    - ``halfplanes`` - a list or tuple of three _MarkedHalfPlanes corresponding
        to the same triangle with vertices z_1, z_2, z_3, with the
        vertices permuted cyclically. I.e. _MarkedHalfPlane(z_1, z_2, z_3), 
        _MarkedHalfPlane(z_2, z_3, z_1) and _MarkedHalfPlane(z_3, z_2, z_1). It 
        would be nicer to have z_1, z_2, z_3 as arguments instead, but much
        slower, because the same _MarkedHalfPlane objects would have to be 
        constructed over and over.

    - ``z`` - a complex number in the unit disk

    OUTPUT:

    - list - a list of three real numbers (the weight-coordinates a, b, c)

    EXAMPLES:

    The following triangle is rotation-symmetric with respect to the
    origin, so the weight-coordianates of the origin are 1/3, 1/3, 1/3:

        >>> z_1 = 0.5
        >>> z_2 = 0.5*complex(cos(2*pi/3), sin(2*pi/3))
        >>> z_3 = z_2.conjugate()
        >>> m_1 = _MarkedHalfPlane(z_1, z_2, z_3)
        >>> m_2 = _MarkedHalfPlane(z_2, z_3, z_1)
        >>> m_3 = _MarkedHalfPlane(z_3, z_1, z_2)
        >>> [round(x, 5) for x in _get_coordinates([m_1, m_2, m_3], 0)]
        [0.33333, 0.33333, 0.33333]

    The coordinates of vertices are 1, 0, 0 in the appropriate order:

        >>> [round(x, 5) for x in _get_coordinates([m_1, m_2, m_3], z_1)]
        [1.0, 0.0, 0.0]
        >>> [round(x, 5) for x in _get_coordinates([m_1, m_2, m_3], z_2)]
        [0.0, 1.0, 0.0]
        >>> [round(x, 5) for x in _get_coordinates([m_1, m_2, m_3], z_3)]
        [0.0, 0.0, 1.0]

    """
    #if z is a vertex, sines_of_angles would lead to division by zero
    for i in range(3):
        if abs(halfplanes[i].z_1 - z) < _epsilon:
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
        if all(weights[i][(i+j)%3] < weights[(i+j)%3][i] + _epsilon for 
                j in {1, 2}):
            lw = i
        
    coord = [weights[lw][i]/weights[i][lw] for i in range(3)]
    return [x/sum(coord) for x in coord]


def _linear_combination(weights, vertices):
    """
    Calculates the linear combination of three pixels with
    given weights.

    INPUT:

    - ``weights`` - list of three real numbers

    - ``vertices`` - list of three tuples of two integers

    OUTPUT:

    - tuple - the tuple of two integers

    EXAMPLES:

        >>> _linear_combination([0.1, 0.2, 0.7], [(0, 0), (0, 100), (100, 0)])
        (70, 20)

    """
    return tuple(_iround(sum(weights[i] * vertices[i][j] for i in range(3))) 
            for j in range(2))

def _complex_to_array_coor(z, halfsize):
    """
    Converts a complex number to array coordinates.

    INPUT:

    - ``z`` - a complex number in the unit disk

    - ``halfsize`` - a float or integer. For an nxn squareshaped image
        it is n/2.

    OUTPUT:

    - tuple of two integers

    EXAMPLES:

        >>> _complex_to_array_coor(0, 500)
        (500, 500)
        >>> _complex_to_array_coor(complex(-0.5, -0.5), 500)
        (750, 250)

    """
    return (_iround(halfsize * (1.0 - z.imag)),_iround(halfsize * (1.0 + z.real)))

def _array_coord_to_complex(pixel, halfsize):
    """
    Converts array coordinates to a complex number.

    INPUT:

    - ``pixel`` - a tuple of two integers

    - ``halfsize`` - a float or integer. For an nxn squareshaped image
        it is n/2.

    OUTPUT:

    - complex number

    EXAMPLES:

        >>> _array_coord_to_complex((0, 0), 500)
        (-1+1j)
        >>> _array_coord_to_complex((750, 250), 500)
        (-0.5-0.5j)

    """
    return complex(pixel[1]/halfsize - 1.0, 1.0 - pixel[0]/halfsize)
    

class TriangleImage(object):
    """
    A triangle-shaped cut from an existing image file.

    TriangleImage objects are used to fill the triangles in the 
    hyperbolic tilings.

    INPUT:

    - ``image_file`` - the name of the image file, as a string. If omitted,
        it is like giving an image file containing one single pixel colored
        to ``default_color`` as input, and setting ``vertices`` to 
        [(0,0), (0,0), (0,0)]. In other words, the hyperbolic triangles
        filled with this TriangleImage will be simply colored to 
        ``default_color``.

    - ``vertices`` - a list of three tuples of two integers (the
        pixel-coordinates of the triangle) or a special string see below). 

        Normally one chooses the 
        coordinates to be valid coordinates for the size of the image,
        but out-of-bounds coordinates are accepted, too. For instance,
        given an image file of size 800x600, the coordinates (0,0),
        (800, 0), (400, 700) bound a triangle where some neighborhood 
        of the vertex (400, 700) is out of the image. Any pixel in the
        triangle that is out-of-bounds is considered to be have color
        ``default_color``.

        Special strings: ``upperleft`` and ``bottomright``. When one of
        these is specified instead of vertices, then the vertice are 
        automatically generated to clip the upper left and bottom right
        triangles of the image.
        

    - ``default_color`` - a tuple of three integers between 0 and 255 as
        the RBG coding of the color. It is white, (255, 255, 255), by
        default.

    EXAMPLES:
        
        >>> t1 = TriangleImage("stinkbug.png", [(0,0), (800, 0), (400, 700)])
        >>> t1 = TriangleImage("stinkbug.png", "upperleft")
        >>> t2 = TriangleImage(default_color = (59, 123, 53))

    """
    def __init__(self, image_file = None, vertices = None, 
            default_color = (255, 255, 255)):
        self.default_color = default_color
        if image_file is None:
            self.array = PIL.Image.new('RGB', (1, 1), default_color).load()
            self.vertices = ((0,0),(0,0),(0,0))
        else:
            im = PIL.Image.open(image_file).convert('RGB')
            self.array = im.load()
            (width, height) = im.size
            topleft = (0, 0)
            topright = (width - 1, 0)
            bottomleft = (0, height - 1)
            bottomright = (width - 1, height - 1)
            if isinstance(vertices, str):
                if vertices == 'upperleft':
                    self.vertices = (topleft, bottomleft, topright)
                elif vertices == 'bottomright':
                    self.vertices = (bottomright, bottomleft, topright)
                else:
                    raise ValueError("Invalid input for ``vertices``.")
            else:
                self.vertices = tuple(tuple(x[::-1]) for x in vertices)
            #self.array = spmisc.imread(image_file).tolist()

def _cosh_side_length(alpha, beta, gamma):
    """
    Apply the law of cosines in hyperbolic geometry to calculate the
    cosh of the length of one side from the angles of the triangle.

    INPUT:

    - ``alpha``, ``beta``, ``gamma`` - the angles of the hyperbolic
        triangle in radians

    OUTPUT:

    - float - the length of the side opposite to ``gamma``

    EXAMPLES:

        >>> _cosh_side_length(pi/8, pi/8, pi/8)
        12.137071184544089

    """    
    return (cos(alpha)*cos(beta) + cos(gamma)) / (sin(alpha) * sin(beta))

def _point_from_origin(cosh_of_dist):
    """
    Returns the positive real number in the Poincare disk if the cosh
    of its hyperbolic distance from the origin is specified.
    
    INPUT:

    - ``cosh_of_dist`` - float, the cosh of the hyperbolic distance
        from the origin

    OUTPUT:

    - float - the positive real number with the property as above

    EXAMPLES:

        >>> _point_from_origin(1.0)
        0.0
        >>> _point_from_origin(2.0)
        0.5773502691896257
        >>> _point_from_origin(3.0)
        0.7071067811865476

    """
    return ((cosh_of_dist - 1.0) / (cosh_of_dist + 1.0))**(0.5)

def _get_triangle(p, q, r):
    """
    Calculates the vertices of a hyperbolic triangle with angles
    pi/p, pi/q, pi/r.

    One of the vertices is the origin, another is on the positive
    real axis, and the third has positive imaginary part.

    INPUT:

    - ``p``, ``q``, ``r`` - positive integers such that 
        1/p + 1/q + 1/r < 1. 

    OUTPUT:

    - tuple of three complex numbers, the vertices of the triangle

    EXAMPLES:

    All the angles of the following triangle are 45 degrees so the
    third vertex is on the line y=x.

        >>> t = _get_triangle(4, 4, 4)
        >>> t[0]
        0.0
        >>> (round(t[2].real, 10), round(t[2].imag, 10))
        (0.4550898606, 0.4550898606)

    """
    if 1.0/p + 1.0/q + 1.0/r >= 1 - _epsilon:
        raise ValueError("The sum of angles in a hyperbolic triangle is less"
                " than pi, so 1/p + 1/q + 1/r < 1 should be satisfied.")
    alpha = pi/p
    beta = pi/q
    gamma = pi/r

    return (0.0, _point_from_origin(_cosh_side_length(alpha, beta, gamma)),
            _point_from_origin(_cosh_side_length(alpha, gamma, beta)) *
            complex(cos(alpha), sin(alpha)))

def create_image(triangle_images, p, q, r, size, output_file_name):
    """
    Writes a new image file with a tiling of the Poincare disk.

    INPUT:

    - ``triangle_images`` - list of TriangleImage objects. Any number of
        TriangleImages may be specified as long as their number divides
        all of p, q, r. (This ensures that the triangles appear
        alternately around each vertex.)

    - ``p``,``q``,``r`` - positive integers such that 
        1/p + 1/q + 1/r < 1/2. They determine how many triangles surround
        the different vertices of each triangle.

    - ``size`` - positive integer. The dimensions of the output file are
        ``size``x``size``. On a 2009 Macbook Pro, a 1000x1000 image takes 
        about 60 seconds to generate with the default CPython interpreter, 
        a little faster (45 seconds) with Cython on PyPy.

    - ``output_file_name`` - string, the name of the output image file. 
        It should be specified with an extension (e.g. "o.jpg", "o.png").
        The allowed formats depend on what libraries are installed with
        PIL.


    EXAMPLES:

    Creating a tiling using one image:

        >>> t = TriangleImage("KillerRabbit.jpg", \
                [(0,0), (800,0), (0,800)]) # doctest: +SKIP
        >>> create_image([t], 4, 8, 10, 1000, "rabbit_circle.jpg") # doctest: +SKIP

    Creating a tiling without a source image, where the triangles alternate
    between red, green, blue:

        >>> t1 = TriangleImage(default_color = (255, 0, 0))
        >>> t2 = TriangleImage(default_color = (0, 255, 0))
        >>> t3 = TriangleImage(default_color = (0, 0, 255))
        >>> create_image([t1, t2, t3], 6, 9, 12, 1000, "colorful.png") # doctest: +SKIP

    """
    n = len(triangle_images)
    vertices = _get_triangle(p, q, r)
    new_image = [[[0]*3 for i in range(size)] for j in range(size)]
    temp_images = [[[[-1]*3 for i in range(size)] for j in range(size)] for x in triangle_images]
    halfplanes = [_MarkedHalfPlane(vertices[i], vertices[(i+1)%3], vertices[(i+2)%3])
            for i in range(3)]
    
    halfsize = size/2.0

    for i in range(size):
        for j in range(size):
            z = _array_coord_to_complex((i, j), halfsize)
            
            if abs(z) > 1 - _epsilon:
                new_image[i][j] = (0, 0, 0)
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
            
            (newi, newj) = _complex_to_array_coor(z, halfsize)
            current_temp_image = temp_images[count % n]
            current_triangle_image = triangle_images[count % n]
            if current_temp_image[newi][newj][0] < 0: #pixel is uninitialized
                weights = _get_coordinates(halfplanes, z)
                (x, y) = _linear_combination(weights, triangle_images[count % n].vertices)
                try:
                    current_temp_image[newi][newj] = current_triangle_image.array[x,y]
                except IndexError: #if the coordinates are out of range
                    current_temp_image[newi][newj] = current_triangle_image.default_color
            new_image[i][j] = current_temp_image[newi][newj]

    #spmisc.imsave(output_file_name, new_image)
    new_img = PIL.Image.new("RGB", (size, size))
    new_img.putdata([x for row in new_image for x in row])
    new_img.save(output_file_name)


def create_image_from_rectangle(input_file_name, p, q, r, size, output_file_name):
    """
    Creates an tiling from a whole rectangular image file.

    In the background it is still a tiling from triangles. There are two
    triangles, the upper left and the lower right triangles of the image.
    Using this function instead of create_image() has the advantage that
    one doesn't have to worry about the vertices of the triangles.

    INPUT:

    - ``input_file_name`` - string, the name of the rectangular image used
        for the tiling.

    - ``p``,``q``,``r`` - positive integers such that 
        1/p + 1/q + 1/r < 1/2. They determine how many triangles surround
        the different vertices of each triangle.

    - ``size`` - positive integer. The dimensions of the output file are
        ``size``x``size``. On a 2009 Macbook Pro, a 1000x1000 image takes 
        about 60 seconds to generate with the default CPython interpreter, 
        a little faster (45 seconds) with Cython on PyPy.

    - ``output_file_name`` - string, the name of the output image file. 
        It should be specified with an extension (e.g. "o.jpg", "o.png").
        The allowed formats depend on what libraries are installed with
        PIL.

    EXAMPLES:

        >>> create_image_from_rectangle("stinkbug.png", 3, 3, 4, \
                100, "output.png") # doctest: +SKIP

    """
    t1 = TriangleImage(input_file_name, 'upperleft')
    t2 = TriangleImage(input_file_name, 'bottomright')
    create_image([t1, t2], p, q, r, size, output_file_name)



if __name__ == "__main__":
    import doctest
    doctest.testmod()

