# Pyreframe, a rudimentary 3D engine made in Python for no real reason.
# By Daniel Rivas, working from various tutorials.

import sys
import pygame
from random import randrange

from engine import (Point, Star, Vector, Camera, LineSegment,
                    Polygon, SymmetricalPolygon, Circle, Renderer, Lighting)

from engine import (scale_space, rotate_space, draw_line,
                    draw_polygon, draw_sympolygon, draw_circle)

# Screen dimensions.
screenWidth = 1000
screenHeight = 1000
desiredDepth = 1000
# Place our light sources.
light0 = Lighting((1000, 0, 0), (255, 0, 0), 1000)
light1 = Lighting((0, 1000, 0), (0, 255, 0), 1000)
light2 = Lighting((0, 0, 1000), (0, 0, 255), 1000)
lights = [light0, light1, light2]

# Define the camera.
camera = Camera(((-1, screenWidth+1),
                 (-1, screenHeight+1),
                 (-1, desiredDepth+1)), lights)
# Build renderer.
renderer = Renderer(camera)


def draw_starfield(renderer):
    # Generate 250 random stars.
    for point in range(250):
        renderer.camera.objectsInWorld.append(
            Star((randrange(0, renderer.screenWidth),
                  randrange(0, renderer.screenHeight),
                  randrange(0, renderer.desiredDepth)), renderer))
    renderer.camera.drawScene()


def draw_pointfield(renderer):
    # Generate 250 random points.
    for point in range(250):
        renderer.camera.objectsInWorld.append(
            Point((randrange(0, renderer.screenWidth),
                   randrange(0, renderer.screenHeight),
                   randrange(0, renderer.desiredDepth)), renderer))
    renderer.camera.drawScene()


def add_content():

    draw_sympolygon(4, Point((500, 500, 100), renderer), 75, renderer)
    draw_sympolygon(4, Point((500, 500, 175), renderer), 75, renderer)
    # draw_circle((500, 500), 50, renderer)

    draw_starfield(renderer)
    draw_pointfield(renderer)


def mainloop(renderer):
    # Main "game" loop.
    print("Press a to scale scene down, s to scale scene up, and r to rotate.")
    add_content()

    while 1:
        renderer.render_screen()
        rotate_space(1,
                     Point((0, 0, 0), renderer),
                     Point((1000, 1000, 1000), renderer), renderer)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_a:
                    scale_space(0.5, renderer)
                elif event.key == pygame.K_s:
                    scale_space(2.0, renderer)
                elif event.key == pygame.K_r:
                    rotate_space(90,
                                 Point((0, 0, 0), renderer),
                                 Point((1000, 1000, 1000), renderer),
                                 renderer)


mainloop(renderer)
