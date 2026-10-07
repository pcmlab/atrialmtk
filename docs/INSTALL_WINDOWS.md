# Installing atrialmtk on Windows

atrialmtk is run on Windows inside **Ubuntu on WSL2** (Windows Subsystem for Linux). Once this guide is complete, every command in the main [README](../README.md) and the example READMEs works unchanged inside the Ubuntu terminal.

Allow about 1–2 hours, most of it waiting for downloads.

**Requirements:** Windows 10 (22H2) or Windows 11, an Intel or AMD (x86-64) processor, at least 8 GB RAM, and about 40 GB of free disk space. Windows on ARM (e.g. Snapdragon "Copilot+ PC" laptops) is not supported, because the openCARP image would run under emulation.

**Two places to type commands.** Some steps are typed in Windows **PowerShell**; most are typed in the **Ubuntu** terminal. Each step says which. Commands starting with `wsl` belong in PowerShell; everything else (`sudo`, `conda`, `git`, `docker`) belongs in Ubuntu.

---

## 1. Prepare Windows

1. Run **Settings → Windows Update** and install all updates, restarting as needed.
2. Check that hardware virtualisation is on: **Task Manager → Performance → CPU** should show *Virtualization: Enabled*. If it says *Disabled*, enable it in your BIOS/UEFI settings (often called *SVM* on AMD or *Intel VT-x*; the key to enter the BIOS at start-up is usually F2, F10 or Del).

## 2. Install WSL2 and Ubuntu

In **PowerShell, run as administrator** (Start → type *PowerShell* → right-click → *Run as administrator*):

```
wsl --install -d Ubuntu-24.04
```

Restart the computer when asked. A *Welcome to Windows Subsystem for Linux* window may appear; you can close it.

Open **Ubuntu 24.04** from the Start menu. The first launch takes a minute, then asks you to create a Linux **username** (lowercase, no spaces) and **password**. Nothing appears while you type the password; this is normal. You will need this password for `sudo` commands.

Check that you are on WSL 2. In **Ubuntu**:

```
uname -r
```

The output should contain `microsoft-standard-WSL2`.

**Opening Ubuntu later:** Start menu → *Ubuntu 24.04*, or type `wsl ~` in PowerShell. (Plain `wsl` starts in the current Windows folder, e.g. `/mnt/c/windows/system32`; type `cd ~` to get to your Linux home folder.)

## 3. Give WSL enough memory

By default WSL uses half of your RAM. On an 8 GB laptop that is too little for Docker plus the pipeline. In **Ubuntu**, run this block as one paste:

```
WINUSER=$(cmd.exe /c "echo %USERNAME%" 2>/dev/null | tr -d '\r')
cat > "/mnt/c/Users/$WINUSER/.wslconfig" <<'CFG'
[wsl2]
memory=6GB
swap=4GB
CFG
cat "/mnt/c/Users/$WINUSER/.wslconfig"
```

(With 16 GB RAM or more, use `memory=10GB` or similar.) If the last line reports *No such file or directory*, run `ls /mnt/c/Users/` to find your Windows profile folder name and use that in place of `$WINUSER`.

Apply it. In **PowerShell**:

```
wsl --shutdown
```

Then reopen Ubuntu and check with `free -h` (the *total* should be close to the value you set).

## 4. Install Docker Desktop

1. Download **Docker Desktop for Windows (AMD64)** from https://www.docker.com/products/docker-desktop/ and run the installer.
2. On the configuration screen keep **Use WSL 2 instead of Hyper-V** ticked. Click OK, then *Close and log out* (or restart) when it finishes.
3. Start **Docker Desktop**, accept the agreement, and choose *Skip* / *Continue without signing in*. Wait until the bottom-left shows **Engine running**.
4. Open **Settings (gear icon) → Resources → WSL integration**, turn on **Ubuntu-24.04**, then click **Apply & restart**.
5. Optional: **Settings → General → Start Docker Desktop when you sign in**.

Close and reopen Ubuntu, then test in **Ubuntu**:

```
docker run --rm hello-world
```

You should see *Hello from Docker!*. If you get `docker: command not found`, check the WSL integration toggle in step 4, run `wsl --shutdown` in PowerShell and reopen Ubuntu.

Docker Desktop must be running whenever you use openCARP. If a docker command hangs or reports *Cannot connect to the Docker daemon*, start Docker Desktop and wait for *Engine running*.

Docker Desktop's licence is free for personal, educational and small-business use; larger organisations may need a subscription. Check with your IT department if unsure.

## 5. Install tools, Miniforge and atrialmtk (in Ubuntu)

All of this section is typed in **Ubuntu**.

```
sudo apt update && sudo apt upgrade -y
sudo apt install -y git build-essential libgl1 libglu1-mesa libxrender1 libxext6
```

Install Miniforge (conda):

```
cd ~
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
exec bash
```

Your prompt should now start with `(base)`.

Clone the repository **into your Linux home folder**:

```
cd ~
git clone https://github.com/pcmlab/atrialmtk.git
cd atrialmtk
```

> Keep atrialmtk and your data inside the Linux home folder (`~`), not under `/mnt/c/...`. Accessing Windows folders from WSL is much slower.

Now create the two conda environments and pull the openCARP image exactly as described in the main README:

- [Conda environments](../README.md#conda-environments): `pointpicking` and `uac`
- [openCARP and meshtool](../README.md#opencarp-and-meshtool): `docker pull docker.opencarp.org/opencarp/opencarp:latest`

## 6. Viewing results

**ParaView (Windows app).** Download the Windows `.msi` from https://www.paraview.org/download/ (the file named `...-Windows-Python3.x-msvc2017-AMD64.msi`, without "MPI") and install it.

To open files that live in WSL, use **File → Open** in ParaView and type this into the file name box, then press Enter:

```
\\wsl$\Ubuntu-24.04\home\<your-linux-username>\atrialmtk
```

You can also browse there in File Explorer and pin the folder to *Quick access*.

**meshalyzer (inside WSL).** Follow the *Linux and Windows (WSL)* meshalyzer instructions in the [main README](../README.md#visualisation-paraview-and-meshalyzer). It runs inside Ubuntu and its window appears on the Windows desktop. If it fails with a FUSE error, use the `--appimage-extract` step described there.

(A native Windows build, `meshalyzer_win64-<version>.zip`, is also available on the meshalyzer releases page. Extract it to a permanent folder such as `C:\Tools\meshalyzer` and run `meshalyzer.exe`. If Windows shows *Windows protected your PC*, click *More info → Run anyway*. This is convenient for viewing results from Windows, but the WSL version is easier to use alongside the pipeline.)

## 7. Check that everything works

In **Ubuntu**:

```
cd ~/atrialmtk
docker run --rm docker.opencarp.org/opencarp/opencarp:latest openCARP -buildinfo | head -3
conda run -n uac python -c "import vtk, numpy, sklearn, meshio; print('uac ok')"
conda run -n pointpicking python -c "import pyvista; print('pyvista ok')"
```

All three should complete without errors. Then follow **[Example-LeftAtrium](../Examples/Example-LeftAtrium/README.md)**.

## Troubleshooting

| Problem | Fix |
|---|---|
| `wsl: command not found` | You typed it inside Ubuntu. `wsl` commands go in PowerShell (or use `wsl.exe` from Ubuntu). |
| Can't find the atrialmtk folder | `cd ~`; plain `wsl` starts in a Windows folder. Use `wsl ~` or open Ubuntu from the Start menu. |
| `docker: command not found` in Ubuntu | Docker Desktop → Settings → Resources → WSL integration → enable Ubuntu-24.04 → Apply & restart; then reopen Ubuntu. |
| `Cannot connect to the Docker daemon` | Start Docker Desktop and wait for *Engine running*. |
| `Unable to locate package libfuse2t64` | Run `sudo apt update` first; on Ubuntu 22.04 the package is `libfuse2`. |
| meshalyzer: FUSE error | Use the `--appimage-extract` method in the main README. |
| Point picking / meshalyzer window doesn't appear | Update WSL in PowerShell with `wsl --update`, then `wsl --shutdown` and reopen Ubuntu. Windows 10 needs build 19044 or later for graphical apps. |
| `conda activate uac` says environment not found | Check `conda env list`. Environment names are case-sensitive on Linux. |
| Everything is very slow | Make sure your files are under `~`, not `/mnt/c/...`, and close other programs while simulations run. |
| `parameter parser error: Unknown parameter` | A newer openCARP release has removed a parameter; see the note in the main README's [openCARP section](../README.md#opencarp-and-meshtool). |
