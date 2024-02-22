---
layout: post
title: "T>T: Compiling DearPyGui on Raspberry Pi 32-bit OS"
date: 2022-06-03
excerpt: "A quick post on how to compile the DearPyGui library on a Raspberry Pi 32-bit OS."
tags: [science, mathematics, programming, GUI, Raspberry, Pi, Raspberry Pi]
comments: false
math: true
---

# The Problem

DearPyGui, available at [https://github.com/hoffstadt/DearPyGui](https://github.com/hoffstadt/DearPyGui), is a robust and speedy Graphical User Interface (GUI) toolkit designed for Python, boasting minimal dependencies. I've been keen to integrate it as the primary library for my AstroPitography GUI. However, a hurdle emerged - it lacks compatibility with the 32-bit Buster Raspberry Pi operating system running on my RPi 3B+.

This blog post serves as a comprehensive guide detailing the steps necessary to compile the DearPyGui library on the 32-bit architecture of a Raspberry Pi, specifically the 3B+. While its primary purpose is to serve as a personal reference for potential future endeavors, I also hope it proves beneficial for others who have found themselves scouring the internet for hours in search of a solution.

# Create a Swap File

Initially I did not worry about swap memory but I ran into `OOM` errors when trying to compile DearPyGui which was a direct consequence of the limited 1GB RAM of the Raspberry Pi 3B+. We are going to create temporary storage space on the SD card using a swap file to act as more RAM for the compilation process. This file is going to swap a section of RAM storage for an idle program and free up memory for other programs.

1. First Create an empty file to use for our swap, **Note** we set the size of the file as 4GiB by: 1K * 4M = 4 GiB.

	```shell
	sudo mkdir -v /var/cache/swap
	cd /var/cache/swap
	sudo dd if=/dev/zero of=swapfile bs=1K count=4M
	sudo chmod 600 swapfile
	```

2. Convert the newly created file into a swap space file.

	```shell
	sudo mkswap swapfile
	```

3. Enable the swap file for paging and swapping.

    ```shell
    sudo swapon swapfile
    ```

4. Add the swap file into the `fstab` file to make it persistent on the next system boot, if required.

    ```shell
    echo "/var/cache/swap/swapfile none swap sw 0 0" | sudo tee -a /etc/fstab
    ```

5. Restart the Raspberry Pi and the new swap file will be ready to use. This can be checked using.

    ```shell
    top -bn1 | grep -i swap
    ```
	and you should see the line.

    ```shell
    KiB Swap:  4194300 total
    ```

# Locally Compile DearPyGui

1. The first step is to set the Raspberry Pi into "Legacy" mode in raspi-config.

2. Install the following libraries using `sudo apt-get install PACKAGE`. Make sure you are also using `python >= 3.6`.

	1. `git`
	2. `cmake`
	3. `libglu1-mesa-dev`
	4. `libgl1-mesa-dev`
	5. `libxrandr-dev`
	6. `libxinerama-dev`
	7. `libxcursor-dev`
	8. `libxi-dev`

3. Clone the DearPyGui repository into a local directory using `git`.

	```shell
	git clone --recursive https://github.com/hoffstadt/DearPyGui
	```

4. Change directory into the new DearPyGui folder directory and run the following.

	```shell
	python3 -m setup bdist_wheel --plat-name linux_armv7l --dist-dir ../dist
	```
	If this is successful a wheel file will be found in `../dist`. **Note** that the platform name, `linux_armv7l` is very important. `pip` considers the name of the resultant wheel file when checking if it can be compiled on the currently available architecture. It will check for `linux_arm7l` when installing.

# Install DearPyGui

We now install the wheel file using pip.

```shell
pip install WHEEL FILE
```

If successful you should now be able to import DearPyGui in Python.

```python
import dearpygui
```