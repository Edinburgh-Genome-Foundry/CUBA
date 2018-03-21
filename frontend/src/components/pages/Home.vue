<!-- src/components/Home.vue -->
<template lang="pug">
.home.page
  h1 The EGF's Collection <br /> of Useful Biological Apps
  web-links(emailSubject='[CUBA] Some user feedback',
            tweetMessage="The Edinburgh Genome Foundry's Collection of Useful Bio Apps",
            tweetUrl='http://cuba.genomefoundry.org')
  .scenarios
    .scenario-category(v-for='category in scenarios', :key='category.category')
      h3 {{category.category}}
      scenariospanel(:scenarios='category.scenarios')

    .egf-codons
      p Powered by
      a(href='https://edinburgh-genome-foundry.github.io/', target='_blank')
        img(src='/static/logos/egf-codons.svg')
</template>

<script>
import scenariospanel from '../widgets/ScenariosPanel'
import scenarios from '../scenarios/scenarios.js'
export default {
  data () {
    return {
      scenarios: scenarios.list
    }
  },
  components: {
    scenariospanel
  }
}
</script>

<style lang='scss' scoped>
h2 {
  margin-bottom: 2em;
  margin-top: 1.5em;
}
p {
  font-size: 18px;
}

.scenarios {
  margin-top: 5em;
}

.egf-codons {
  width: 100%;
  margin-top: 10em;
  text-align: center;
  img {
    margin: 0 auto;
    width: 140px;
  }
}
</style>
