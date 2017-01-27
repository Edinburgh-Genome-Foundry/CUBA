<template lang="pug">
el-row
  el-col(:xs="24", :sm="12", :md="6", :lg="6" v-for='scenario in scenarios')
    el-card.box-card
      img(:src='scenario.infos.icon', alt='')
      .caption
        h3 {{scenario.infos.title}}
        p.thumbnail-description {{scenario.infos.description}}
        router-link.btn.btn-primary(:to='scenario.infos.path') Go
</template>

<script>
import scenarios from './scenarios'
export default {
  data: () => ({
    scenarios: scenarios.list
  })
}
</script>

<style scoped>


.box-card {
  width:90%;
  margin: 0px auto;
  text-align: center
}
img {
  height: 80px
}

.thumbnail .caption {
  padding: 0px;
}

.thumbnail-description {
  font-size: 14px;
}

</style>
