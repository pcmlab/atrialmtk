# atrialmtk
**Atrial Modelling Toolkit for Constructing Bilayer and Volumetric Atrial Models at Scale**

To enable large in silico trials and personalised model predictions on clinical timescales, it is imperative that models can be constructed quickly and reproducibly. First, we aimed to overcome the challenges of constructing cardiac models at scale through developing a robust, open-source pipeline for bilayer and volumetric atrial models. Second, we aimed to investigate the effects of fibres, fibrosis and model representation (surface or volumetric) on fibrillatory dynamics. To construct bilayer and volumetric models, we extended our previously developed coordinate system to incorporate transmurality, atrial regions, and fibres (rule-based or data driven diffusion tensor MRI). We demonstrate the generalisation of these methods by creating a cohort of 1000 biatrial bilayer and volumetric models derived from CT data, as well as models from MRI, and electroanatomic mapping. Fibrillatory dynamics diverged between bilayer and volumetric simulations across the CT cohort (correlation co-efficient for phase singularity maps: LA 0.27±0.19, RA 0.41±0.14). Adding fibrotic remodelling stabilised re-entries and reduced the impact of model type in the LA (LA: 0.52±0.20, RA: 0.36±0.18). The choice of fibre field has a small effect on paced activation data (<12ms), but a larger effect on fibrillatory dynamics. Overall, we developed an open-source user-friendly pipeline for generating atrial models from imaging or electroanatomic mapping enabling in silico clinical trials at scale.

![Model pipeline](https://github.com/pcmlab/atrialmtk/blob/main/images/Figure1_Schematicv2.jpg?raw=true)

# Installation instructions

atrialmtk runs on **Linux**, **macOS** (Intel and Apple Silicon) and **Windows** (via WSL2). You will need:

| Component | Used for | Install |
|---|---|---|
| Docker | Runs openCARP and meshtool | [below](#docker) |
| openCARP + meshtool | Laplace solves and simulations | [below](#opencarp-and-meshtool) |
| Miniforge (conda) | Python environments `pointpicking` and `uac` | [below](#conda-environments) |
| ParaView and/or meshalyzer | Viewing meshes and results | [below](#visualisation-paraview-and-meshalyzer) |

> **Windows users:** follow the step-by-step guide in **[docs/INSTALL_WINDOWS.md](docs/INSTALL_WINDOWS.md)**. Everything runs inside Ubuntu on WSL2, and the commands in the rest of this README then work unchanged.

## Docker

- **macOS / Windows:** install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and start it. Docker Desktop must be running whenever you use openCARP.
- **Linux:** install [Docker Engine](https://docs.docker.com/engine/install/), then allow your user to run docker without `sudo`:

    ```
    sudo groupadd docker
    sudo usermod -aG docker $USER
    newgrp docker
    ```

    Log out and back in (or restart) for this to take effect.

Check that docker works:

```
docker run --rm hello-world
```

## openCARP and meshtool

Pull the openCARP image (this includes meshtool):

```
docker pull docker.opencarp.org/opencarp/opencarp:latest
docker run --rm docker.opencarp.org/opencarp/opencarp:latest openCARP -buildinfo
```

The second command should print the openCARP version. To install openCARP directly (for example on a HPC system) see https://opencarp.org/download/installation.

> **Note on openCARP versions:** recent openCARP releases removed the parameters `ellip_use_pt`, `parab_use_pt` and `mat_entries_per_row`. These are commented out in the `.par` files in this repository. If you see `parameter parser error: Unknown parameter ...` with a future openCARP release, comment out the named parameter in the `.par` file and please open an issue.
>
> **Apple Silicon Macs:** the openCARP image runs under emulation, which works but is slower. Use short simulations when testing.

## Conda environments

We recommend [Miniforge](https://github.com/conda-forge/miniforge), which uses the conda-forge channel and a fast solver. 

Install Miniforge (Linux, WSL and macOS):

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
~/miniforge3/bin/conda init "$(basename "$SHELL")"
```

Then open a new terminal. The environments only need to be created once.

**1a. Point picking environment (`pointpicking`):**

```
conda create -n pointpicking -c conda-forge python=3.10 pandas numpy -y
conda activate pointpicking
python -m pip install pyvista==0.42.2 vtk==9.2.6
conda deactivate
```

Run the point picking code with `python Rough_Point_Picking.py` while this environment is active.

**1b. Universal Atrial Coordinates environment (`uac`):**

Create a Python 3.8 environment with conda, then install the pinned packages with pip (this takes about a minute; solving the full environment with conda can take a very long time). From the top level of this repository:

```
conda create -n uac -c conda-forge python=3.8 -y
conda activate uac
pip install -r src/3Processing/UAC_Codes/requirements.txt
```

**Apple Silicon Macs (M1 and later)** must create this environment as an Intel (x86) environment, because `vtk==9.0.3` has no native arm64 build:

```
CONDA_SUBDIR=osx-64 conda create -n uac -c conda-forge python=3.8 -y
conda activate uac
conda config --env --set subdir osx-64
pip install -r src/3Processing/UAC_Codes/requirements.txt
```

Check the environment:

```
conda activate uac
python -c "import vtk, numpy, sklearn, meshio; print('uac ok', vtk.VTK_VERSION, numpy.__version__)"
```

This should print `uac ok 9.0.3 1.23.1`. With openCARP installed you can now run the processing scripts (e.g. `./mri-la.sh` from `src/3Processing`). Type `conda deactivate` when finished.

> If you created this environment previously under the name `UAC`, either keep using `conda activate UAC` or remove it with `conda env remove -n UAC` and recreate it as above.

## Visualisation: ParaView and meshalyzer

**ParaView** (recommended on all platforms): download from https://www.paraview.org/download/

- Windows: the `.msi` installer named `...-Windows-Python3.x-msvc2017-AMD64.msi` (without "MPI" in the name).
- macOS: the `.dmg` for your chip (`arm64` for Apple Silicon, `x86_64` for Intel). Check via Apple menu, About This Mac.
- Linux: the `.tar.gz`; extract it and run `bin/paraview`.

**meshalyzer**: releases are at https://git.opencarp.org/openCARP/meshalyzer/-/releases. Builds are provided for Linux and Windows only; on macOS it must be built from source.

*Linux and Windows (WSL):*

```
sudo apt update
sudo apt install -y libfuse2t64 libgl1 libglu1-mesa || sudo apt install -y libfuse2 libgl1 libglu1-mesa
mkdir -p ~/bin
mv ~/Downloads/Meshalyzer-*-x86_64.AppImage ~/bin/meshalyzer    # adjust to where you saved it
chmod +x ~/bin/meshalyzer
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
meshalyzer
```

If it fails with a FUSE error (common on WSL), extract the AppImage instead:

```
cd ~/bin
./meshalyzer --appimage-extract && mv squashfs-root meshalyzer-app && rm meshalyzer
ln -s ~/bin/meshalyzer-app/AppRun ~/bin/meshalyzer
```

*macOS (build from source):*

```
brew install cmake fltk glew freetype libpng pkg-config libomp
git clone https://git.opencarp.org/openCARP/meshalyzer.git
cd meshalyzer
conda deactivate          # repeat until no environment is active, so conda libraries are not picked up
make -j
ln -sf "$PWD/meshalyzer" /opt/homebrew/bin/meshalyzer    # Intel Macs: /usr/local/bin/meshalyzer
```

Build from a git clone rather than the source zip: the build needs the `.git` folder.

*Updating meshalyzer:*

- Linux / WSL: download the new AppImage and replace `~/bin/meshalyzer` with it (or, if you used the extracted version, delete `~/bin/meshalyzer-app` and repeat the extract step).
- macOS: `cd` to your meshalyzer clone, then `git pull && make clean && make -j`. The symlink picks up the new build automatically.
- Check the version with `meshalyzer --version` (or the About window).

# **Usage** 

We have included the following examples: 

## **[Example-LeftAtrium](/Examples/Example-LeftAtrium/README.md)**
This example generates a left atrial surface model from the Atrial Challenge dataset.  All steps are explained in the following categories:

0ImagingData

1Clipping

2Landmarks

3Processing

4Simulation

![LA model](https://github.com/pcmlab/atrialmtk/blob/main/images/laP.png?raw=true)


## **[Example-Biatrial-MRI](/Examples/Example-Biatrial-MRI/README.md)**

This example generates a biatrial bilayer model from the Atrial Challenge dataset. This example can be run after the Example-LeftAtrium to add the right atrial component of the model, and interatrial connections as a biatrial bilayer model. The steps are detailed in the same categories as for Example-LeftAtrium. 

![Biatrial model](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial3.png?raw=true)

## **[Example-Biatrial-CT-shape-model](/Examples/Example-Biatrial-CT-shape-model/README.md)**

This example generates a biatrial bilayer model and a biatrial volumetric model from a CT-derived statistical shape model. In this case, we start with the atrial surfaces, and progress through the steps to generate two versions of the model in a bilayer format (triangular mesh) and a volumetric format (tetrahedral mesh). The steps are as follows:

0SurfaceMeshData 

2Landmarks

3Processing

4Simulation

(these meshes do not require any clipping). 

![Biatrial CT model](https://github.com/pcmlab/atrialmtk/blob/main/images/ct-biatrial.png?raw=true)
