<template lang='pug'>
div
  div(v-for='(selector, index) in value')
    el-row
      el-col(:span=22)
        digestionselector(v-model='value[index]')
      el-col(:span=2)
        el-button(@click='value.splice(index, 1)', icon='delete')
  .addbutton
    el-button(@click='value.push([])' icon='plus') Add a digestion
</template>

<script>
import digestionselector from './DigestionSelector'
export default {
  props: {
    value: {
      default: () => ([[]])
    }
  },
  watch: {
    value: function (val) {
      this.$emit('input', val)
    }
  },
  components: { digestionselector }

}
</script>

<style scoped>
.addbutton {
  text-align: center;
}
</style>
