Escher
======

Python script is for creating Escher-like tilings of the Poincare disk.

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
