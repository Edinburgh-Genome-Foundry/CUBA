<template lang="pug">
el-select(:value='value', placeholder='Enzymes in this digestion', filterable, :multiple='true')
  el-option(v-for='enzyme in enzymes.list', :label='enzyme', :value='enzyme')
</template>

<script>
import enzymes from '../../assets/enzymes'

export default {
  props: {
    value: {
      default: () => ([])
    }
  },
  data: function () {
    return {
      enzymes: enzymes
    }
  },
  watch: {
    value: function (val) {
      this.$emit('input', val)
    }
  }
}
</script>

<style scoped>
.el-tag {
  height: 35px !important;
}

.el-select {
  width: 100%;
}
</style>
