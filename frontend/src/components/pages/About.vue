<template lang='pug'>
.page
  h1  About this website
  img.about-big-logo(src='../../assets/images/cuba-title.png')
  p
    |  The Edinburgh Genome Foundry offers DNA assembly and other lab operations as
    |  a service, and develops software to support DNA design and assembly.
    |  This website offers some functionalities which we thought may be of interest
    |  to other foundries or bench scientists. Do not hesitate to&nbsp
    a(href="mailto:egf-software@ed.ac.uk?subject=[CUBA] Comment or request&body=Hey EGF,") contact us
    | &nbspfor comments or requests.
</template>

<style scoped>
.about-big-logo {
  width: 50%;
  margin-left: 25%;
}
p {
  text-align: justify;
  width: 60%;
  max-width: 600px;
  margin: 0 auto;
}
</style>
