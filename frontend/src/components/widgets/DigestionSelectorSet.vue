<template lang='pug'>
div
  div(v-for='(selector, index) in value')
    el-row(:gutter='20')
      el-col(:span=22)
        digestionselector(v-model='value[index]')
      el-col(:span=2)
        el-button(@click='value.splice(index, 1)', icon='el-icon-delete' circle)
  .addbutton(style='margin-top: 2em')
    el-button(@click='value.push([])' icon='el-icon-plus') Add a digestion
</template>

<script>
import digestionselector from './DigestionSelector'
export default {
  props: {
    value: {
      default: () => ([[]])
    }
  },
  data () {
    return {
      digestions: this.value
    }
  },
  watch: {
    digestions: {
      deep: true,
      handler (val) {
        this.$emit('input', val)
      }
    },
    value: {
      deep: true,
      handler (val) {
        this.digestions = val
      }
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
