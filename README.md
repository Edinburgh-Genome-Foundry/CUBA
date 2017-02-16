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

Installing and deploying
-------------------------

First you need to install ``docker`` and ``docker-compose`` on your machine. That's the difficult part.

Then install download your CAB application like this:

```
git clone CAB.git
```

Go to the root folder (the one containing this README and the ``docker-compose.yml``) and type

```
docker-compose build
```

This may take some time, especially if it is the first time you build a CAB application on your machine.

To launch your application in development mode, just type

```
docker-compose up
```

Now you can access the app from your browser on port 8080 of your machine, either
[http://localhost:8080](http://localhost:8080/) if it is running on your local
machine or something like  ``http://139.59.162.184:8080`` if the CAB app is running
on a distant server.

At this stage, any change you do to your code, back-end or front-end, will be reflected live without you needed to refresh.


To run in production mode (this will turn down debug mode in the frontend and backend).

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up
```
