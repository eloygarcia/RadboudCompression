#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 16:42:05 2021

@author: eloy
"""
import os.path
from typing import Optional
import docker
from subprocess import call
from shutil import copyfile


#from compDistances import computeDistances

CONTAINER_NAME = 'NIFTYSIM-CONTAINER'
tempFolder = '/home/eloygarcia/Escritorio/Pruebas/temporalNiftySim'

# https://dev.to/serhatteker/how-to-check-if-a-docker-container-running-with-python-3aoj
def is_container_running(container_name: str) -> Optional[bool]:
    """Verify the status of a container by it's name

    :param container_name: the name of the container
    :return: boolean or None
    """
    RUNNING = "running"
    # Connect to Docker using the default socket or the configuration
    # in your environment
    docker_client = docker.from_env()
    # Or give configuration
    # docker_socket = "unix://var/run/docker.sock"
    # docker_client = docker.DockerClient(docker_socket)

    ### ELOY
    ## check if niftysim image exists
    images = docker_client.images.list(filters = { "reference" : "niftysim:latest"})
    if len(images)==0:
        print("Does not exist latest version of the niftysim image!")
        return False

    ## chek if container exists and it is running
    containers = docker_client.containers.list(all=True, filters={'name':CONTAINER_NAME})
    if len(containers) == 0:
        container = docker_client.containers.run("niftysim:latest", name=CONTAINER_NAME,
                                                 volumes={tempFolder: {'bind': '/home/results', 'mode': 'rw'}},
                                                 detach=True, tty=True, stdin_open=False, runtime='nvidia')
        print('Container created successfully!!')

    try:
        container = docker_client.containers.get(CONTAINER_NAME)
    except docker.errors.NotFound as exc:
        print(f"Check container name!\n{exc.explanation}")
        return False
    else:
        container_state = container.attrs["State"]

    if not container_state["Status"] == RUNNING:
        print('Container is not running')
        container.restart()
    return container

    """
    try:
        container = docker_client.containers.get(container_name)
    except docker.errors.NotFound as exc:
        print(f"Check container name!\n{exc.explanation}")
        return False
    else:
        container_state = container.attrs["State"]
        return container_state["Status"] == RUNNING
    """

def compressionFunction(patient, thickness, mesh_path, output_path, gravity):
    ### Create NiftySim Model
    text = ['/home/eloygarcia/Escritorio/Pruebas/Release Compression/FEM Simulation/Write Nifty Model/WriteNiftyModel',
            mesh_path, output_path, str(thickness), str(gravity)]
    print(text)
    call(text)

    ### Exec Niftysim into docker container
    copyfile(output_path, os.path.join(tempFolder, os.path.basename(output_path)))
    container = is_container_running(CONTAINER_NAME)
    print('is here')

    text = ['home/niftysim/release/source/niftysim',
        '-x', os.path.join('/home/results', os.path.basename(output_path)),
        '-v', '-t', '-sport',
        '-export-mesh',
         os.path.join('/home/results', patient+'-outputMesh.vtk')]
    print('Beggining NiftySim compression')
    container.exec_run(' '.join(text))
    print('Finishing compression')

    print( os.path.join(tempFolder, patient+'-outputMesh.vtk') )
    print( os.path.join(os.path.dirname(output_path), patient+'-outputMesh.vtk') )
    copyfile(os.path.join(tempFolder, patient+'-outputMesh.vtk'),
             os.path.join(os.path.dirname(output_path), patient+'-outputMesh.vtk'))

    ### Reconstructing Compressed mesh

    text = ['/home/eloygarcia/Escritorio/Pruebas/Release Compression/Phantom Reconstruction/Extract Compressed Mesh/getCompressedMesh',
            os.path.join(os.path.dirname(output_path), patient+'-outputMesh.vtk'),
            os.path.join(os.path.dirname(output_path), patient+'-'+str(gravity)+'-compressedMesh.vtk')]
    call(text)