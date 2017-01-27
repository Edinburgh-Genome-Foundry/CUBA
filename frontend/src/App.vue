<template lang='pug'>
  #app
    .container
      navbar
      router-view
    footersection
</template>



<script>
import auth from './auth'
import navbar from './components/Navbar'
import footersection from './components/Footer'

export default {
  data () {
    return {
      user: auth.user
    }
  },
  methods: {
    logout () {
      auth.logout()
    }
  },
  components: {
    navbar,
    footersection
  }
}
</script>

<style>
#app {
  font-family: 'Raleway', 'Source Sans Pro',  Helvetica, Arial, sans-serif;
  font-size:16px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  color: #2c3e50;
  margin-top: 60px;
}

html {
  position: relative;
  min-height: 100%;
}

body {
  /* Margin bottom by footer height */
  margin-bottom: 60px;
}

.container {
  width:850px;
  max-width:95%;
  margin: 0px auto;
}



h1, h2 {
  text-align: center;
  margin-bottom:1.5em

}




.router-link-active {
    color: "black";
}
</style>
