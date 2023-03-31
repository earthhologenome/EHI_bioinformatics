# EHI_Bioinformatics 
# 🐨->💩->🦠->🧬->🖥️->😏
Bioinformatics pipeline to process EHI data.

*updated 29/03/2023, Raphael Eisenhofer*

#### General information:
This pipeline uses [![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io), and manages dependencies using conda (or mamba) for reproducibility and deployability. The 0_Code directory contains the snakefiles, scripts, and conda environment yamls. 

#### Getting started:
Firstly, you'll need to set up an alias for connecting to ERDA -- this is **essential** for proper function.

If you haven't already, create a public ssh key:
```
ssh-keygen
```

This will prompt you with something like this:
```
Generating public/private rsa key pair.
Enter file in which to save the key (/home/ncl550/.ssh/id_rsa):
```

Press enter.

Next, we need to setup your ssh config file. If this doesn't exist already, create it here (**using your own user name, of course!**):
```
/home/ncl550/.ssh/config
```

And within this file, add your erda alias like so (**use your email address**):
```
Host erda ucph-erda
 Hostname io.erda.dk
 VerifyHostKeyDNS yes
 User raphael.eisenhofer@sund.ku.dk
 Port 22
```

Now copy your public key (you'll need to change your path to match your user number):
```
more /home/ncl550/.ssh/id_rsa.pub
```

Copy the large string (e.g.)
```
ssh-rsa asd79as7d98as7d987as987d98a79d87a897d9a7d98a98d9 ncl550@mjolnirhead01fl.unicph.domain
```

Now, log into ERDA via your browser `https://erda.dk/wsgi-bin/home.py`

Navigate to your 'settings' in the bottom right corner:

![setup1](figures/setup1.png)

And click the SFTP tab at the top:

![setup1](figures/setup2.png)

You then want to paste your public SSH key into the box (see yellow arrow). Then click save.

Now, relog into mjolnir (or reload your shell if you know how to), and test the alias:
```
sftp erda
```

It may prompt you to allow the connection, if so, press enter.

Great, you can now launch the preprocessing pipeline via the 'PR_script' column in the 'PR BATCH' tab on AirTable (**Launch view**)! Note that you should select your email address in the 'Email' column to get alerted upon completion. 

## WARNING
### ERDA only allows 16 concurrent sessions per user, so for now, please only launch a maximum of 2 preprocessing batches at any given time.

The pipeline has been written such that the output stats are automatically inserted into the EHI AirTable through the magical powers of API, so sit back and relax.

Voila!
