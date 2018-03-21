
<p align="center">
<img alt="DNA Cauldron Logo" title="DNA Cauldron Logo" src="https://github.com/Edinburgh-Genome-Foundry/CUBA/raw/master/frontend/src/assets/images/cuba-title.png" width="400">
</p>
<h1 align="center"> The EGF Collection of Useful Biological Apps </h1>

[![Build Status](https://travis-ci.org/Edinburgh-Genome-Foundry/CUBA.svg?branch=master)](https://travis-ci.org/Edinburgh-Genome-Foundry/CUBA)

This repository contains the source code of [cuba.genomefoundry.org](http://cuba.genomefoundry.org/),
a website of the Edinburgh Genome Foundry enabling anyone to use some of the EGF's
biological software.

![screenshots](https://github.com/Edinburgh-Genome-Foundry/CUBA/raw/master/docs/imgs/screenshots.png)




## How is it built ?

CUBA is based on the [CAB](https://github.com/Edinburgh-Genome-Foundry/CAB)
boilerplate, making it easy to create new apps with a form in the frontend and
custom computations in the backend. It features job
queues (with progress feedback for the user), form widgets like file uploaders,
help buttons, and many more.

### Getting started

The next steps will download, install, and launch CUBA on your computer.

1. Install ``docker`` and ``docker-compose`` on your machine. This step depends
   on your machine (Windows, Linux, MacOS) so you'll need to google it.

2. Download CUBA from Github:

```
git clone git+github.com/Edinburgh-Genome-Foundry/CUBA.git
```

3. Go to the root ``CUBA/`` directory (the one containing this README and the
   ``docker-compose.yml``) and launch your application in development mode. The
   first time you try this, Docker will download and build a lot of things,
   which may take several minutes. It will only take a few seconds the next
   times you run this command.

```
docker-compose up
```

4. Go to your browser and type ``localhost`` or ``127.0.0.1`` in the address bar.
   You should see the website appear. the console in which you launched
   ``docker-compose`` will keep printing logs of the different components
   (django, vue) so you can keep track and debug.

### Creating a new app

The next steps will add a new app to the CUBA project.

1. Go to ``frontend/src/components/scenarios`` and create a new scenario view
   with a form, for instance by duplicating the file
   [ExampleScenario.vue](https://github.com/Edinburgh-Genome-Foundry/CUBA/blob/master/frontend/src/components/scenarios/ExampleScenario.vue).

2. Register your scenario in file ``scenarios.js`` (in the same folder)
   by adding ``require('./ExampleScenario')`` under the category you want.
   You should now see your new scenario in the home page and the menu of the
   website.

3. Next we will add some backend computations to process the form and return a
   result. First go to ``backend/app/views`` and create a new folder
   on the model of [``/example_scenario``](https://github.com/Edinburgh-Genome-Foundry/CUBA/tree/master/backend/app/views/example_scenario).

4. Register the scenario in ``backend/app/views/__init__.py`` by adding

```
from .example_scenario import ExampleScenarioView
```

5. Register the URL by adding the following line at the end of
   ``backend/website/urls.py``:

```
url(r'^api/start/example_scenario$',
      views.ExampleScenarioView.as_view()),
```

6. That's it. You now have a new app with frontend and backend !

### Deploying the website on the web

The next steps will put your website on the web. Note that many other deployment
workflows are possible.

1. Get a hosting server (for instance from Amazon Web Services or Digital Ocean).
   Get the IP address of this server (for instance ``123.12.123.123`` but
   it could be ``mydomain.com`` if you have registered this domain and pointed it
   to your server).

2. Log into this server (``ssh root@123.12.123.123`` or ``ssh root@mydomain.com``)
   and install Docker and Docker-Compose (some Digital Ocean servers come with
   these already installed).

3. From your computer, in the CUBA root directory, run the following command to
   create a code repository on the distant server, and register that distant
   repository under the name ``prodserver``

```
./init_remote_git.sh root@123.12.123.123 CUBA prodserver
```

4. On the remote server, in the folder ``CUBA.git``, start the website in
   production mode:

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up --build -d
```

5. Wait some time and go in your browser at the address ``123.12.123.123``,
   your website should be live !

Every time you want to update the website, from your computer in the CUBA root
directory run ``git push prodserver master``. You will need to rebuild the
containers on the server if you have modified the frontend or added dependencies
to the backend (we may simplify this later).

Licence
-------

CUBA is an open source software originally written by [Zulko](https://github.com/Zulko)
at the [Edinburgh Genome Foundry](http://genomefoundry.org/) and released on
[Github](https://github.com/Edinburgh-Genome-Foundry/CUBA) under the MIT licence
(¢ Edinburgh Genome Foundry). Everyone is welcome to contribute !


More biology software ?
-----------------------
<p>
<a href="https://edinburgh-genome-foundry.github.io/">
   <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png" />
</a>
</p>

EGF Codon is powered by the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_
synthetic biology software suite for DNA design, manufacturing and validation.
