#!/usr/bin/env python3
"""
Undirected Graph Visualizer using Pygame

Reads graph from stdin:ute.html
Not enough edge pairs supplied (expected 250 edges).
First line: N M
Next M lines: u v

Node numbering can be 0-based (0..N-1) or 1-based (1..N). The program auto-detects and adjusts.

Controls:
 - Drag nodes with left mouse button
 - Space: pause / resume force layout
 - r: reset random positions
 - q or ESC or window close: quit

Run:
    pip install pygame
    python graph_visualizer.py < input.txt

"""
import sys
import math
import random
import pygame
from collections import defaultdict

# ---------------------
# Input parsing
# ---------------------


def read_graph():
    data = sys.stdin.read().strip().split()
    if not data:
        print("No input provided. Expecting: N M then M pairs u v", file=sys.stderr)
        sys.exit(1)
    it = iter(data)
    try:
        N = int(next(it))
        M = int(next(it))
    except StopIteration:
        print("Input too short. Expecting N M", file=sys.stderr)
        sys.exit(1)

    edges = []
    nodes_seen = set()
    for _ in range(M):
        try:
            u = int(next(it))
            v = int(next(it))
        except StopIteration:
            print("Not enough edge pairs supplied (expected {} edges).".format(
                M), file=sys.stderr)
            sys.exit(1)
        edges.append((u, v))
        nodes_seen.add(u)
        nodes_seen.add(v)

    # detect 1-based numbering if any index >= N
    if nodes_seen and max(nodes_seen) >= N:
        # assume 1-based input, convert to 0-based
        edges = [(u-1, v-1) for (u, v) in edges]

    # clamp and validate
    for (u, v) in edges:
        if not (0 <= u < N and 0 <= v < N):
            print(f"Edge ({u},{v}) out of range for N={N}", file=sys.stderr)
            sys.exit(1)

    return N, edges

# ---------------------
# Force-directed layout
# ---------------------


class GraphLayout:
    def __init__(self, N, edges, width, height):
        self.N = N
        self.edges = edges
        self.width = width
        self.height = height

        self.pos = [pygame.math.Vector2(random.uniform(
            50, width-50), random.uniform(50, height-50)) for _ in range(N)]
        self.vel = [pygame.math.Vector2(0, 0) for _ in range(N)]

        self.adj = [[] for _ in range(N)]
        for u, v in edges:
            if v not in self.adj[u]:
                self.adj[u].append(v)
            if u not in self.adj[v]:
                self.adj[v].append(u)

        # parameters
        self.k_attract = 0.01
        self.k_repulse = 10000.0
        self.damping = 0.85
        self.max_disp = 10.0

    def step(self):
        # compute repulsive forces
        forces = [pygame.math.Vector2(0, 0) for _ in range(self.N)]
        for i in range(self.N):
            for j in range(i+1, self.N):
                delta = self.pos[i] - self.pos[j]
                dist_sq = delta.x*delta.x + delta.y*delta.y + 1e-6
                dist = math.sqrt(dist_sq)
                # repulsive magnitude
                f = self.k_repulse / dist_sq
                direction = delta / dist
                forces[i] += direction * f
                forces[j] -= direction * f

        # attractive forces along edges
        for u, v in self.edges:
            delta = self.pos[v] - self.pos[u]
            dist = max(delta.length(), 1e-6)
            # spring-like attraction towards target length
            # desired length depends on window size
            desired = 120.0
            diff = dist - desired
            f = self.k_attract * diff
            direction = delta / dist
            forces[u] += direction * f
            forces[v] -= direction * f

        # integrate
        for i in range(self.N):
            acc = forces[i]
            self.vel[i] = (self.vel[i] + acc) * self.damping
            # limit displacement to avoid blowups
            if self.vel[i].length() > self.max_disp:
                self.vel[i].scale_to_length(self.max_disp)
            self.pos[i] += self.vel[i]

            # keep inside window with soft bounds
            margin = 40
            if self.pos[i].x < margin:
                self.pos[i].x = margin + (self.pos[i].x - margin) * 0.5
                self.vel[i].x *= -0.3
            if self.pos[i].x > self.width - margin:
                self.pos[i].x = (self.width - margin) + \
                    (self.pos[i].x - (self.width - margin)) * 0.5
                self.vel[i].x *= -0.3
            if self.pos[i].y < margin:
                self.pos[i].y = margin + (self.pos[i].y - margin) * 0.5
                self.vel[i].y *= -0.3
            if self.pos[i].y > self.height - margin:
                self.pos[i].y = (self.height - margin) + \
                    (self.pos[i].y - (self.height - margin)) * 0.5
                self.vel[i].y *= -0.3

    def reset_positions(self):
        self.pos = [pygame.math.Vector2(random.uniform(
            50, self.width-50), random.uniform(50, self.height-50)) for _ in range(self.N)]
        self.vel = [pygame.math.Vector2(0, 0) for _ in range(self.N)]

# ---------------------
# Pygame drawing and interaction
# ---------------------


def run_visualizer(N, edges):
    pygame.init()
    WIDTH, HEIGHT = 1920, 1080
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    pygame.display.set_caption('Undirected Graph Visualizer')
    clock = pygame.time.Clock()

    layout = GraphLayout(N, edges, WIDTH, HEIGHT)

    node_radius = 12
    font = pygame.font.SysFont(None, 18)

    running = True
    paused = False
    dragging = None
    drag_offset = pygame.math.Vector2(0, 0)

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    paused = not paused
                elif event.key == pygame.K_r:
                    layout.reset_positions()
                elif event.key == pygame.K_q or event.key == pygame.K_ESCAPE:
                    running = False
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:  # left click -> pick node
                    mx, my = event.pos
                    best = None
                    best_dist = 1e9
                    for i, pos in enumerate(layout.pos):
                        d = (pos.x - mx)**2 + (pos.y - my)**2
                        if d < best_dist and d <= (node_radius*node_radius*4):
                            best = i
                            best_dist = d
                    if best is not None:
                        dragging = best
                        drag_offset = layout.pos[best] - \
                            pygame.math.Vector2(mx, my)
            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:
                    dragging = None
            elif event.type == pygame.MOUSEMOTION:
                if dragging is not None:
                    mx, my = event.pos
                    layout.pos[dragging] = pygame.math.Vector2(
                        mx, my) + drag_offset
                    layout.vel[dragging] = pygame.math.Vector2(0, 0)

        if not paused:
            # perform several small steps per frame for stability
            for _ in range(2):
                layout.step()

        screen.fill((30, 30, 30))

        # draw edges
        for u, v in edges:
            pygame.draw.aaline(screen, (180, 180, 180),
                               layout.pos[u], layout.pos[v])

        # draw nodes
        for i, pos in enumerate(layout.pos):
            # circle fill
            color = (100, 190, 255) if i == dragging else (120, 200, 150)
            pygame.draw.circle(
                screen, color, (int(pos.x), int(pos.y)), node_radius)
            # border
            pygame.draw.circle(screen, (20, 20, 20),
                               (int(pos.x), int(pos.y)), node_radius, 2)
            # label
            label = font.render(str(i), True, (250, 250, 250))
            label_rect = label.get_rect(center=(int(pos.x), int(pos.y)))
            screen.blit(label, label_rect)

        # UI overlay
        help_lines = [
            "Left-click + drag: move node",
            "Space: pause/resume layout",
            "r: randomize positions",
            "q or ESC: quit",
            f"Nodes: {N}    Edges: {len(edges)}",
        ]
        for idx, line in enumerate(help_lines):
            surf = font.render(line, True, (220, 220, 220))
            screen.blit(surf, (10, 10 + idx*18))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()

# ---------------------
# main
# ---------------------


def main():
    N, edges = read_graph()
    run_visualizer(N, edges)


if __name__ == '__main__':
    main()
