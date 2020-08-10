---
layout: post
title: "T>T: Langton's Ant 2: Electric Boogaloo"
date: 2020-07-26
excerpt: "Further investigation into the complex emergent behaviour of langton's ant."
tags: [science, mathematics, programming, matplotlib, script, langtons, ant, python]
comments: false
math: true
---

In a previous post we programmed a simple implementation of Langton's ant which followed the movement pattern:

1) If you are on a white square, turn 90° clockwise, flip the colour of the square, then move forward one square.

2) If you are on a black square, turn 90° anti-clockwise, flip the colour of the square, then move forward one square.

This movement pattern is classed as the `RL` scheme, as the first step involves turning right, whilst the second step involves turning left, alternating the colour each time. There are an infinity of possible extensions to this idea; and we will explore some of these in this post including:
Greater than binary colour scheme
1) More than two colour options for steps.

2) Multiple ants on the same grid.

We will create a generalised program for any movement pattern, e.g. `LRRL` and for any number of ants; and then have fun testing out different combinations.

# `RL` Langton's Ant
