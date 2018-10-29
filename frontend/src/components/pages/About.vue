<template lang='pug'>
.page
  img.about-big-logo(src='../../assets/images/cuba-title.png')
  h4(align='center') About us
  p.
    The Edinburgh Genome Foundry offers DNA assembly and other lab operations as
    a service, and develops software to support DNA design and assembly.
    This website offers some functionalities which we thought may be of interest
    to other foundries or bench scientists. Do not hesitate to&nbsp
    #[a(href="mailto:egf-software@ed.ac.uk?subject=[CUBA] Comment or request&body=Hey EGF,") contact us]
    &nbspfor comments or requests.
  
  h4(style='text-align: center; margin-top: 3em') Privacy notice
  p.
    We do not record DNA sequences, DNA schemas, or any other form data
    sent to our servers. This data is only kept on our servers for the time
    it takes to complete the computations, and up to 3 minutes after the
    result is sent to you.

  p.
    We do record metadata associated with requests (such as request url,
    user IP address, or connection time) in order to detect improper website
    usage (bots etc.) 
</template>

<style scoped>
.about-big-logo {
  width: 50%;
  margin-left: 25%;
  margin-top: 3em;
}
p {
  text-align: justify;
  width: 60%;
  max-width: 600px;
  margin: 1em auto;
}
</style>
