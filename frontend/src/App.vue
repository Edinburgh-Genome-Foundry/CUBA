<template lang='pug'>
  #website
    github-corners(repo="Edinburgh-Genome-Foundry" bg-color='#8ac1e8')
    .page
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

<style lang='scss'>

@import './assets/css/animate.css';
@import url(https://fonts.googleapis.com/css?family=Open+Sans);
@import url(https://fonts.googleapis.com/css?family=Inconsolata);
@import url(https://fonts.googleapis.com/css?family=Raleway);

#website {
  font-family: 'Raleway', 'Source Sans Pro',  Helvetica, Arial, sans-serif;
  font-size:16px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  color: #2c3e50;
  margin-top: 60px;
}
.el-select-dropdown, .el-input__inner {
  font-family: 'Raleway', 'Source Sans Pro',  Helvetica, Arial, sans-serif;
}
html {
  position: relative;
  min-height: 100%;
}

body {
  /* Margin bottom by footer height */
  margin-bottom: 60px;
}

.page, .navbar {
  width:850px;
  max-width:95%;
  margin: 0px auto;
}

hr {
    border : 0;
    height: 1px;
    margin-top: 60px;
    margin-bottom: 60px;
    background-image: linear-gradient(to right, rgba(0, 0, 0, 0), rgba(0, 0, 0, 0.75), rgba(0, 0, 0, 0));
}

.el-button {
  font-size: 18px;
}

.el-button.center {
  display: block;
  margin: 2em auto;
}

h1, h2 {
  text-align: center;
  margin-bottom:0.5em;
  font-weight: normal
}

h1 { font-size: 3em; }
h2 { font-size: 2em; }

img.icon {
  display: block;
  margin: 0 auto;
  max-width: 100px;
}

.center {
  text-align: center;
}


.router-link-active {
    color: "black";
}

.scenario-description {
  max-width: 400px;
  margin: 1em auto 1em;
}

.el-popper {
  word-break: normal !important;
}
.el-steps {
  max-width:400px;
  margin: 0 auto;
  background: none !important;
}
.el-steps--simple{
  background: none !important;
}

p.inline {
  span.options{
    margin-left: 1em;
  }
}

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 60px;
}

.el-checkbox.inline {
  margin-left: 15px;
}


.el-slider.inline {
  display: inline-block;
  margin-left: 20px;
  margin-right: 20px;
  margin-bottom: -12px;
}

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

.el-tooltip__popper {
  font-family: 'Raleway' !important;
}
</style>
