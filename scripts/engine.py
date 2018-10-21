# Pyreframe, a rudimentary 3D engine made in Python for no real reason.
# By Daniel Rivas.

import math
import pygame


class Point(object):
    def __init__(self, coords, renderer):
        """
        Define point's XYZ coordinates as three-tuple.
        self.num = coords
        self.renderer = renderer
        self.colour = 255, 255, 255
        """

    def addVectorToPoint(self, vals):
        # Add a vector to this point, resulting in a point.
        num1, num2, num3 = self.num
        val1, val2, val3 = vals.num

        num1 += val1
        num2 += val2
        num3 += val3

        return Point((num1, num2, num3), self.renderer)

    def subtractVectorFromPoint(self, vals):
        """
        Subtract a vector from this point, resulting in a point.
        """

        num1, num2, num3 = self.num
        val1, val2, val3 = vals.num

        num1 -= val1
        num2 -= val2
        num3 -= val3

        return Point((num1, num2, num3), self.renderer)

    def subtractPointFromPoint(self, vals):
        """
        Construct a vector from this point to another point.
        """
        num1, num2, num3 = self.num
        val1, val2, val3 = vals.num

        num1 -= val1
        num2 -= val2
        num3 -= val3

        return Vector((num1, num2, num3), self.renderer)

    def drawPoint(self, colour):
        """
        Use Renderer's function to draw this.
        """
        self.renderer.render_point(self.num, colour)

    def setPointToPoint(self, newpos):
        self.num = newpos.num


class Star(Point):
    """
    Like the Point() class, but rendered using circles to indicate depth.
    Makes for a nice visual effect.
    """

    def drawPoint(self, colour):
        """
        Draw using main rendering technique.
        """
        self.renderer.render_star(self.num, colour)


class Vector(object):
    def __init__(self, coords, renderer):
        """
        Define this vector's XYZ position with a three-tuple.
        """
        self.num = (coords)
        self.renderer = renderer

    def addVectorToVector(self, vals):
        """
        Add a vector to this vector, resulting in a vector.
        """
        num1, num2, num3 = self.num
        val1, val2, val3 = vals.num

        num1 += val1
        num2 += val2
        num3 += val3

        return Vector((num1, num2, num3), self.renderer)

    def subtractVectorFromVector(self, vals):
        """
        Subtract a vector from this vector, resulting in a vector.
        """
        num1, num2, num3 = self.num
        val1, val2, val3 = vals.num

        num1 -= val1
        num2 -= val2
        num3 -= val3

        return Vector((num1, num2, num3), self.renderer)

    def rotateVectorXY(self, deg):
        """
        Rotate this vector along the XY axis.
        Convert degrees input to radians.
        """
        rad = math.radians(deg)
        # Rotation matrix for XY, split into 3 tuples for convenience.
        transform1 = (math.cos(rad), -math.sin(rad), 0)
        transform2 = (math.sin(rad), math.cos(rad), 0)
        transform3 = (0, 0, 1)
        # Transform.
        num1, num2, num3 = self.num
        out1 = ((num1 * transform1[0]) + (num2 * transform1[1]) +
                (num3 * transform1[2]))

        out2 = ((num1 * transform2[0]) + (num2 * transform2[1]) +
                (num3 * transform2[2]))

        out3 = ((num1 * transform3[0]) + (num2 * transform3[1]) +
                (num3 * transform3[2]))

        return Vector((out1, out2, out3), self.renderer)

    def rotateVectorXZ(self, deg):
        """
        Rotate this vector by the XZ axis.
        Convert degrees input to radians.
        """
        rad = math.radians(deg)
        # Rotation matrix for XZ, split into 3 tuples for convenience.
        transform1 = (math.cos(rad), 0, math.sin(rad))
        transform2 = (0, 1, 0)
        transform3 = (-math.sin(rad), 0, math.cos(rad))
        # Transform
        num1, num2, num3 = self.num
        out1 = ((num1 * transform1[0]) + (num2 * transform1[1]) +
                (num3 * transform1[2]))

        out2 = ((num1 * transform2[0]) + (num2 * transform2[1]) +
                (num3 * transform2[2]))

        out3 = ((num1 * transform3[0]) + (num2 * transform3[1]) +
                (num3 * transform3[2]))

        return Vector((out1, out2, out3), self.renderer)

    def rotateVectorYZ(self, deg):
        """
        Rotate this vector by the YZ axis.
        Convert degrees input to radians.
        """
        rad = math.radians(deg)
        # Rotation matrix for YZ, split into 3 tuples for convenience.
        transform1 = (1, 0, 0)
        transform2 = (0, math.cos(rad), -math.sin(rad))
        transform3 = (0, math.sin(rad), math.cos(rad))
        # Transform
        num1, num2, num3 = self.num
        out1 = ((num1 * transform1[0]) + (num2 * transform1[1]) +
                (num3 * transform1[2]))

        out2 = ((num1 * transform2[0]) + (num2 * transform2[1]) +
                (num3 * transform2[2]))

        out3 = ((num1 * transform3[0]) + (num2 * transform3[1]) +
                (num3 * transform3[2]))

        return Vector((out1, out2, out3), self.renderer)

    def rotateVectorByAxis(self, p1, p2, theta):
        """
        Rotates a vector around an arbitrary axis.
        Translate the space's origin to the axis's origin.
        """
        p = (self.num[0] - p1.num[0], self.num[1] - p1.num[1],
             self.num[2] - p1.num[2])
        # Initialise point q, which we will rotate.
        q = Point((0, 0, 0), self.renderer)
        # Vector matrix N describes the axis we rotate around.
        N = (p2.num[0] - p1.num[0], p2.num[1] - p1.num[1],
             p2.num[2] - p1.num[2])
        # Scalar value Nm is the length of vector N.
        Nm = math.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
        # Unit vector n.
        n = Vector((N[0] / Nm, N[1] / Nm, N[2] / Nm), self.renderer)

        # Matrix common factors, for readability of rotation matrix.
        c = math.cos(theta)
        t = (1 - math.cos(theta))
        s = math.sin(theta)
        X = n.num[0]
        Y = n.num[1]
        Z = n.num[2]

        # Rotation matrix.
        d11 = t * X**2 + c
        d12 = t * X * Y - s * Z
        d13 = t * X * Z + s * Y
        d21 = t * X * Y + s * Z
        d22 = t * Y**2 + c
        d23 = t * Y * Z - s * X
        d31 = t * X * Z - s * Y
        d32 = t * Y * Z + s * X
        d33 = t * Z**2 + c

        # Multiply position matrix p by rotation matrix M,
        # and set q's coordinates to the result.
        q.setPointToPoint(
            Point((d11 * p[0] + d12 * p[1] + d13 * p[2],
                   d21 * p[0] + d22 * p[1] + d23 * p[2],
                   d31 * p[0] + d32 * p[1] + d33 * p[2]), self.renderer))
        # Translate axis and rotated point back to original location.
        q.setPointToPoint(
            Point((q.num[0] + p1.num[0], q.num[1] + p1.num[1],
                   q.num[2] + p1.num[2]), self.renderer))
        return q

    def scaleVector(self, scalePoint, scale):
        """
        Scales the vector by an input three-tuple, around an arbitrary point.
        """
        num0, num1, num2 = self.num
        # Translate space so the scalePoint lies on the origin.
        p0 = num0 - scalePoint.num[0]
        p1 = num1 - scalePoint.num[1]
        p2 = num2 - scalePoint.num[2]
        # Scale the point relative to the origin.
        p0 *= scale[0]
        p1 *= scale[1]
        p2 *= scale[2]
        # Move space back again.
        return Vector((p0 + scalePoint.num[0], p1 + scalePoint.num[1],
                       p2 + scalePoint.num[2]), self.renderer)


class Camera(object):
    def __init__(self, bounds, lights):
        """
        Define minimum and maximum bounds.
        Bounds should be a tuple containing three two-tuples.
        """
        self.minX, self.maxX = bounds[0]
        self.minY, self.maxY = bounds[1]
        self.minZ, self.maxZ = bounds[2]
        # An array containing all existing objects.
        self.objectsInWorld = []
        # Array containing Lighting objects, our light sources.
        self.lights = lights

    def drawScene(self):
        """
        Draw every object that is within the camera bounds.
        """
        objcounter = 0
        # Draw each point that is within bounds.
        for obj in self.objectsInWorld:
            if (type(obj) is Point or type(obj) is Star):
                if (obj.num[0] > self.minX and obj.num[0] < self.maxX
                        and obj.num[1] > self.minY and obj.num[1] < self.maxY
                        and obj.num[2] > self.minZ and obj.num[2] < self.maxZ):
                    objcounter += 1
                    # Check what colour the point should be,
                    # given the light sources present in the scene.
                    colour = self.checkPointColour(obj)
                    # Don't bother drawing invisible points.
                    if colour != [0, 0, 0]:
                        # Draw the point.
                        obj.drawPoint(colour)
            elif (type(obj) is LineSegment or type(obj) is Circle):

                segment_array = obj.returnPoints()
                for point in segment_array:
                    if (point.num[0] > self.minX and point.num[0] < self.maxX
                            and point.num[1] > self.minY
                            and point.num[1] < self.maxY
                            and point.num[2] > self.minZ
                            and point.num[2] < self.maxZ):
                        colour = self.checkPointColour(point)
                        point.drawPoint(colour)
            elif (type(obj) is Polygon or type(obj) is SymmetricalPolygon):

                new_lines = obj.returnLines()
                self.objectsInWorld += new_lines

    def checkPointColour(self, point):
        """
        Takes a point and decides what colour it should be,
        Based on that point's proximity to each light source in the scene.
        """
        red = 0
        green = 0
        blue = 0
        light_num = 0
        for light in self.lights:
            light_num += 1

            # Calculate distance between point and the light
            d_v = [
                light.num[0] - point.num[0], light.num[1] - point.num[1],
                light.num[2] - point.num[2]
            ]
            # Pythagoras's theorem.
            distance = math.sqrt(d_v[0]**2 + d_v[1]**2 + d_v[2]**2)

            # Check if light illuminates the point at all.
            if distance <= light.radius:
                # Then calculate the illumination percentage and apply it.
                percentage = distance / light.radius
                red += (light.colour[0] * percentage)
                green += (light.colour[1] * percentage)
                blue += (light.colour[2] * percentage)

        colour = red, green, blue
        return colour


class Lighting(object):
    """
    This class represents a light source.
    Place in the scene and it will illuminate any points within its
    radius, based on the point's distance from the source.
    """

    def __init__(self, num, colour, radius):
        self.num = num
        self.colour = colour
        self.radius = radius


class LineSegment(object):
    def __init__(self, start, end, renderer):
        """
        Define starting and ending points from two 2-tuples.
        """

        self.start = start
        self.end = end
        self.renderer = renderer

    def returnPoints(self):
        """
        Uses Bresenham's algorithm (3D version) to find the points
        between our starting and ending points.
        """
        # Empty array for holding our line's points.
        pointArray = []
        # Assign start/end points to variables for ease of use.
        x0, y0, z0 = self.start.num
        x1, y1, z1 = self.end.num

        x0 = int(x0)
        x1 = int(x1)
        y0 = int(y0)
        y1 = int(y1)
        z0 = int(z0)
        z1 = int(z1)

        # Vector differences dx, dy and dz.
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        dz = abs(z1 - z0)
        # Step by 1 pixel each time.
        sx = math.copysign(1, x1 - x0)
        sy = math.copysign(1, y1 - y0)
        sz = math.copysign(1, z1 - z0)
        # Add the first point.
        pointArray.append(Point((x0, y0, z0), self.renderer))

        # The driving axis is the one with the greatest distance between
        # starting and ending points.

        # If x is the driving axis:
        if (dx >= dy and dx >= dz):
            err1 = 2 * dy - dx
            err2 = 2 * dz - dx
            while x0 != x1:
                x0 += sx
                if err1 >= 0:
                    y0 += sy
                    err1 -= 2 * dx
                if err2 >= 0:
                    z0 += sz
                    err2 -= 2 * dx
                err1 += 2 * dy
                err2 += 2 * dz
                pointArray.append(Point((x0, y0, z0), self.renderer))

        # If y is the driving axis:
        if (dy >= dx and dy >= dz):
            err1 = 2 * dx - dy
            err2 = 2 * dz - dy
            while y0 != y1:
                y0 += sy
                if err1 >= 0:
                    x0 += sx
                    err1 -= 2 * dy
                if err2 >= 0:
                    z0 += sz
                    err2 -= 2 * dy
                err1 += 2 * dx
                err2 += 2 * dz
                pointArray.append(Point((x0, y0, z0), self.renderer))

        # If z is the driving axis:
        else:
            err1 = 2 * dy - dz
            err2 = 2 * dx - dz
            while z0 != z1:
                z0 += sz
                if err1 >= 0:
                    y0 += sy
                    err1 -= 2 * dz
                if err2 >= 0:
                    x0 += sx
                    err2 -= 2 * dz
                err1 += 2 * dy
                err2 += 2 * dx
                pointArray.append(Point((x0, y0, z0), self.renderer))

        return pointArray


class Polygon(object):
    def __init__(self, points, renderer):

        self.points = points
        self.renderer = renderer

    def returnLines(self):

        returning_lines = []
        for point in range(len(self.points)):
            if point == len(self.points) - 1:
                returning_lines.append(
                    LineSegment(self.points[point], self.points[0],
                                self.renderer))
            else:
                returning_lines.append(
                    LineSegment(self.points[point], self.points[point + 1],
                                self.renderer))

        return returning_lines


class SymmetricalPolygon(Polygon):
    def __init__(self, vertices, centre, radius, renderer):
        self.vertices = vertices
        self.centre = centre
        self.radius = radius
        self.renderer = renderer
        self.points = self.findPoints()

    def findPoints(self):
        point_array = []
        tempVector = Vector((0, self.radius, 0), self.renderer)
        angle = 360 / self.vertices

        for vertex in range(self.vertices):
            point_array.append(self.centre.addVectorToPoint(tempVector))
            tempVector = tempVector.rotateVectorXY(angle)

        return point_array


class Circle(object):
    """
    Draws a circle of arbitrary size and location.
    """

    def __init__(self, centre, radius, renderer):
        self.centre = centre
        self.radius = radius
        self.renderer = renderer

    def returnPoints(self):
        """
        Returns an array of points describing a circle with given centre
        and radius.
        """
        pointArray = []
        # Values needed for the algorithm.
        # f tracks the progress.
        f = 1 - self.radius
        xstep = 1
        ystep = -2 * self.radius
        x = 0
        y = self.radius

        # Algorithm doesn't build maxima, so we do it manually.
        pointArray.append(
            Point((self.centre[0], self.centre[1] + self.radius, 0),
                  self.renderer))
        pointArray.append(
            Point((self.centre[0], self.centre[1] - self.radius, 0),
                  self.renderer))
        pointArray.append(
            Point((self.centre[0] + self.radius, self.centre[1], 0),
                  self.renderer))
        pointArray.append(
            Point((self.centre[0] - self.radius, self.centre[1], 0),
                  self.renderer))
        # Build the rest of the points algorithmically.
        while x < y:
            if f >= 0:
                y -= 1
                ystep += 2
                f += ystep
            x += 1
            xstep += 2
            f += xstep

            # Build the current arc.
            pointArray.append(
                Point((self.centre[0] + x, self.centre[1] + y, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] - x, self.centre[1] + y, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] + x, self.centre[1] - y, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] - x, self.centre[1] - y, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] + y, self.centre[1] + x, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] - y, self.centre[1] + x, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] + y, self.centre[1] - x, 0),
                      self.renderer))
            pointArray.append(
                Point((self.centre[0] - y, self.centre[1] - x, 0),
                      self.renderer))

        return pointArray


class Renderer(object):
    def __init__(self, camera):
        """
        Initialise Pygame.
        Screen/3D-space dimensions.
        """
        self.camera = camera
        self.screenWidth = self.camera.maxX
        self.screenHeight = self.camera.maxY
        self.desiredDepth = self.camera.maxZ

        pygame.init()
        self.size = [self.screenWidth, self.screenHeight]
        self.black = 0, 0, 0
        self.colour = 255, 255, 255

        self.screen = pygame.display.set_mode(self.size)

    def render_point(self, location, colour):
        """
        Draw a pixel at the Point's coordinates X/Y.
        """
        colour = (int(colour[0]), int(colour[1]), int(colour[2]))
        self.screen.set_at((int(location[0]), int(location[1])),
                           (colour))  # <-- Colour.

    def render_star(self, location, colour):
        """
        Draw a circle onscreen to represent a point.
        Z-position represented by circle-size.
        """
        pointsize = 4 - int(10 * (location[2] / self.desiredDepth))
        # Check to see if the point is in front or behind.
        # If it's behind the viewer, don't draw it.
        # (Also pygame crashes if given a negative circle radius.)
        colour = (int(colour[0]), int(colour[1]), int(colour[2]))
        if pointsize > -1:
            pygame.draw.circle(self.screen, colour,
                               (int(location[0]), int(location[1])), pointsize)

    def render_screen(self):
        pygame.display.flip()
        self.screen.fill(self.colour)


def scale_space(amount, renderer):
    """
    Scale every point by the same amount, relative to the origin.
    """
    tempVector = Vector((0, 0, 0), renderer)
    origin = Point((0, 0, 0), renderer)
    # scalePoint is the (arbitrary) point we're scaling relative to.
    # In this case, the centre of the simulated space.
    scalePoint = Point((renderer.screenWidth / 2, renderer.screenHeight / 2,
                        renderer.desiredDepth / 2), renderer)

    for obj in renderer.camera.objectsInWorld:
        if type(obj) is Point or type(obj) is Star:
            tempVector = obj.subtractPointFromPoint(origin)
            tempVector = tempVector.scaleVector(scalePoint,
                                                (amount, amount, amount))
            obj.setPointToPoint(tempVector)
        elif type(obj) is LineSegment:
            tempVector1 = obj.start.subtractPointFromPoint(origin)
            tempVector1 = tempVector1.scaleVector(scalePoint,
                                                  (amount, amount, amount))

            tempVector2 = obj.end.subtractPointFromPoint(origin)
            tempVector2 = tempVector2.scaleVector(scalePoint,
                                                  (amount, amount, amount))

            obj.start.setPointToPoint(tempVector1)
            obj.end.setPointToPoint(tempVector2)

    renderer.camera.drawScene()


def rotate_space(degree, p1, p2, renderer):
    """
    Rotate every point by the same amount about an axis.
    """
    tempVector = Vector((0, 0, 0), renderer)
    tempVector1 = Vector((0, 0, 0), renderer)
    tempVector2 = Vector((0, 0, 0), renderer)
    origin = Point((0, 0, 0), renderer)
    theta = math.radians(degree)

    for obj in renderer.camera.objectsInWorld:
        if type(obj) is Point or type(obj) is Star:
            tempVector = obj.subtractPointFromPoint(origin)
            tempVector = tempVector.rotateVectorByAxis(p1, p2, theta)
            obj.setPointToPoint(tempVector)
        elif type(obj) is LineSegment:
            tempVector1 = obj.start.subtractPointFromPoint(origin)
            tempVector1 = tempVector1.rotateVectorByAxis(p1, p2, theta)
            tempVector2 = obj.end.subtractPointFromPoint(origin)
            tempVector2 = tempVector2.rotateVectorByAxis(p1, p2, theta)
            obj.start.setPointToPoint(tempVector1)
            obj.end.setPointToPoint(tempVector2)

    renderer.camera.drawScene()


def draw_line(start, end, renderer):
    """
    Generate a line of points between two Points.
    """
    myLine = LineSegment(start, end, renderer)
    renderer.camera.objectsInWorld.append(myLine)
    renderer.camera.drawScene()


def draw_polygon(points, renderer):
    """
    Generate a polygon from arbitrary number of Points.
    """
    myPolygon = Polygon(points, renderer)
    renderer.camera.objectsInWorld += myPolygon.returnLines()
    renderer.camera.drawScene()


def draw_sympolygon(vertices, centre, radius, renderer):
    """
    Auto-generate a symmetrical polygon.
    """
    myPolygon = SymmetricalPolygon(vertices, centre, radius, renderer)
    renderer.camera.objectsInWorld += myPolygon.returnLines()
    renderer.camera.drawScene()


def draw_circle(location, radius, renderer):
    myCircle = Circle(location, radius, renderer)
    renderer.camera.objectsInWorld.append(myCircle)
    renderer.camera.drawScene()
