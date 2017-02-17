CAB - The friendly Computational Applications Boilerplate
----------------------------------------------------------

CAB makes it easy to do webpages for the following scenario: the user connects
to the website and lands on a page listing all the apps (or "scenarios"). The user selects an app,
select some options and parameters, click on "send". Then the computer thinks for a while
(while giving progress feedback) then return a report which can be HTML (pictures, tables) or/and a file (PDF, ZIP...)

CAB was written to

Why use CAB ?
-------------

CAB only uses tools that are both powerful and fun to work with:

- Python and Django in the backend, to run any scientific/computational software.
- Vue.js in the frontend, one of the easiest frameworks to learn and maintain.
- Docker-compose for deploying your app anywhere in just a few lines.

This comes with the following advantages:

- Complete live-reload ! The changes you make to your backend or frontend code
  take effect immediately. No need to refresh a page or restart a server by hand.
- Plenty of easy-to-use library of components to build your user interface without headache
  (by default, CAB uses [Element](http://element.eleme.io/#/en-US))
- Will work on any system, one-line install, one-line deploy, and without any
  clash with your current settings as everything runs inside a container.

## Installing and developing

First you need to install ``docker`` and ``docker-compose`` on your machine. That's the difficult part.

Then install download your CAB application like this:

```
git clone CAB.git
```

Go to the root folder (the one containing this README and the ``docker-compose.yml``) and type this to launch your application in development mode

```
docker-compose up
```

Note that the first time you type this, Docker will need to build the necessary
containers, which may take a few minutes. Then this command just takes second.

Now you can access the app from your browser on port 8080 of your machine, either
[http://localhost:8080](http://localhost:8080/) if it is running on your local
machine or something like  ``http://139.59.162.184:8080`` if the CAB app is running
on a distant server.

At this stage, any change you do to your code, back-end or front-end, will be reflected live without you needed to refresh.

## Use in production


To run in production mode (this will turn down debug mode in the frontend and backend).

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```




```
git remote add production git+ssh:123.456.789.123:
```

```
git push production master
```

Licence
-------

CAB is an open source software originally written by [Zulko](https://github.com/Zulko)
at the [Edinburgh Genome Foundry](http://genomefoundry.org/) and released on [Github](#) under the MIT licence (¢ Edinburgh Genome Foundry). Everyone is welcome to contribute.

If you publish an app made with CAB, you can licence your own code under any other licence and copyright,  by placing the terms of the licence on top of each new file you create (this will typically be one file in the frontend and one file in the backend).
