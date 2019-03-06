<template lang="pug">
el-select(v-model='ladder', placeholder='Select')
  el-option(v-for='item in ladder_options', :label='item.label',
            :value='item.value', :key='item.value')
</template>

<script>
export default {
  props: {
    value: {default: '100_to_10k'}
  },
  data () {
    return {
      ladder: this.value,
      ladder_options: [
        {
          label: 'Ladder 100 bp - 4000 bp',
          value: '100_to_4k'
        },
        {
          label: 'Ladder 35 bp - 5 kbp (AATI)',
          value: '35_to_5k'
        },
        {
          label: 'Ladder 75 bp - 15 kbp (AATI)',
          value: '75_to_15k'
        },
        {
          label: 'Ladder 100 bp - 10 kbp (NEB 2-log)',
          value: '100_to_10k'
        },
        {
          label: 'Ladder 500 bp - 10 kbp (NEB 1k)',
          value: '500_to_10k'
        },
        {
          label: 'Ladder 100 bp - 1350 bp (NEB 50bp)',
          value: '100_to_1350'
        }
      ]
    }
  },
  watch: {
    ladder () {
      this.$emit('input', this.ladder)
    }
  }
}
</script>

<style scoped>
</style>
